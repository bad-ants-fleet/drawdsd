#
# High-level and mid-level library interface.
#
import logging
log = logging.getLogger(__name__)
import math
import numpy as np

from dsdobjects.objectio import set_io_objects, clear_io_objects, read_pil, read_pil_line
from itertools import combinations

from .rendering import get_drawing, get_rgb_palette, draw_stem, draw_tentacles, draw_arrowheads
from .components import fourway_module, hairpin_module

def agl(a):
    return a % 360

def almost_same(point1, point2):
    if math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2) < 1:
        return True
    else:
        return False

def get_default_plot_params(cplx):
    """ Produce the defaults for:
        - pair_angles
        - loop_lengths
    """
    ptable = list(cplx.pair_table)
    stable = list(cplx.strand_table)
    pair_angles = []
    loop_lengths = []
    for si, strand in enumerate(ptable):
        nlen, llen = [], 0
        for di, pair in enumerate(strand):
            if pair is None: # unpaired
                llen += stable[si][di].length
            else: # pair
                if pair > (si, di):
                    pair_angles.append(0)
                nlen.append(llen)
                llen = 0
        nlen.append(llen)
        loop_lengths.append(nlen)
    return pair_angles, loop_lengths

def get_segments(ptable, loop_lengths):
    """Helper function to decompose a pair table into segments.

    Segments[(pair, pair)] = ([segment1, segment2], (ll1, ll2))

    """
    segments = {}
    identity, seg, pos = None, [[], []], 0
    ll1, ll2, lli = 0, 0, 0 # loop_length stuff
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pos:
                ll2 = loop_lengths[si][lli]
            else:
                ll1 = loop_lengths[si][lli]
            if pair is not None:
                pos += 1
                lli += 1
                if identity is not None:
                    segments[identity] = seg, (ll1, ll2)
                    identity, seg, pos = None, [[], []], 1
                    ll1, ll2 = 0, 0
                identity = (si, di), pair
            else:
                seg[pos].append([si, di])
        # Strand break means we start a new segment!
        segments[identity] = seg, (ll1, ll2)
        identity, seg, pos = None, [[], []], 0
        ll1, ll2, lli = 0, 0, 0
    return segments

def get_raw_modules(ptable, stable, segments, pair_angles, scale):
    modules, metadata = [], []
    # We start at 5' end, so the modules must be in order!
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pair is None or (si, di) > pair:
                continue
            (sj, dj) = pair
            # Get both segments by their identifier.
            sg1, (k1, k2) = segments[(si, di), (sj, dj)]
            sg2, (k3, k4) = segments[(sj, dj), (si, di)]
            # The stem is formed via the single paired domain.
            stem = [stable[si][di], stable[sj][dj]] 
            # The tentacles listed in segments.
            t1 = [stable[i][j] for (i, j) in sg1[0]]
            t4 = [stable[i][j] for (i, j) in sg2[1]]
            # Check if it is a hairpin!
            if sg1[1] and sg2[0] == [] and sg1[1][-1] == [sj, dj-1]:
                th = [stable[i][j] for (i, j) in sg1[1]]
                m = hairpin_module(stem, t1, t4, th, scale = scale)
                # Provide the loop_length data
                if k1: m.k1 = scale * k1
                if k4: m.k4 = scale * k4
                # We pretend it has four ends to save some code later.
                meta = [[*sg1[0], [si, di], *sg1[1]], 
                        [[sj, dj], *sg2[1]]]
            else:
                t2 = [stable[i][j] for (i, j) in sg1[1]]
                t3 = [stable[i][j] for (i, j) in sg2[0]]
                m = fourway_module(stem, t1, t2, t3, t4, scale = scale)
                # Provide the loop_length data
                if k1: m.k1 = scale * k1
                if k2: m.k2 = scale * k2
                if k3: m.k3 = scale * k3
                if k4: m.k4 = scale * k4
                meta = [[*sg1[0], [si, di], *sg1[1]], 
                        [*sg2[0], [sj, dj], *sg2[1]]]
            # Provide the pair_angle data
            m.angle = pair_angles.pop(0)
            modules.append(m)
            metadata.append(meta)
    return modules, metadata

def get_final_modules(cplx, pair_angles, loop_lengths, origin = (0, 0), scale = 10):
    stable = list(cplx.strand_table)
    ptable = list(cplx.pair_table)

    # Get the helper datastructure segments.
    segments = get_segments(ptable, loop_lengths)

    # Returns all modules but has no knowledge on their neighborhood.
    modules, metadata = get_raw_modules(ptable, stable, segments, pair_angles, scale)

    # Now set all the small angles
    def helper_21(m, n):
        assert m._nma2 is None
        assert n._pma1 is None
        m._nma2 = n.angle
        n._pma1 = m.angle
        assert n.mp1[0] is n
        assert m.mp2[0] is m
        n.mp1 = (m, 'p2')
        m.mp2 = (n, 'p1')

    def helper_41(m, n):
        assert m._nma4 is None
        assert n._pma1 is None
        m._nma4 = n.angle
        n._pma1 = agl(m.angle + 180)
        assert n.mp1[0] is n
        assert m.mp4[0] is m
        n.mp1 = (m, 'p4')
        m.mp4 = (n, 'p1')

    def helper_43(m, n):
        assert m._nma4 is None
        assert n._pma3 is None
        m._nma4 = agl(n.angle + 180)
        n._pma3 = agl(m.angle + 180)
        assert m.mp4[0] is m
        assert n.mp3[0] is n
        m.mp4 = (n, 'p3')
        n.mp3 = (m, 'p4')

    for (m, md), (n, nd) in combinations(zip(modules, metadata), 2):
        # Get all tentacle ends.
        mt1, mt2 = md[0][0], md[0][-1] # first and last
        mt3, mt4 = md[1][0], md[1][-1]
        nt1, nt2 = nd[0][0], nd[0][-1]
        nt3, nt4 = nd[1][0], nd[1][-1]
        # The strand continues from m(2) to n(1)
        if (mt2[0], mt2[1]) == (nt1[0], nt1[1] - 1):
            helper_21(m, n)
        # The strand continues from n(2) to m(1)
        if (nt2[0], nt2[1]) == (mt1[0], mt1[1] - 1):
            helper_21(n, m)
        # The strand continues from m(4) to n(1)
        if (mt4[0], mt4[1]) == (nt1[0], nt1[1] - 1):
            helper_41(m, n)
        # The strand continues from n(4) to m(1)
        if (nt4[0], nt4[1]) == (mt1[0], mt1[1] - 1):
            helper_41(n, m)
        # The strand continues from m(4) to n(3)
        if (mt4[0], mt4[1]) == (nt3[0], nt3[1] - 1):
            helper_43(m, n)
        # The strand continues from n(4) to m(3)
        if (nt4[0], nt4[1]) == (mt3[0], mt3[1] - 1):
            helper_43(n, m)
        # The strand continues from n(2) to m(3)
        if (nt2[0], nt2[1]) == (mt3[0], mt3[1] - 1):
            raise NotImplementedError('Wait, do we have a pseudoknot?')
        # The strand continues from m(2) to n(3)
        if (mt2[0], mt2[1]) == (nt3[0], nt3[1]-1):
            raise NotImplementedError('Wait, do we have a pseudoknot?')

    # Finally, infer all points.
    m = modules[0]
    m.p1 = origin
    for m in modules:
        m.angle = m.angle # Re-setting the angle updates all small angles!
        m.p1 = getattr(*m.mp1)
        if hasattr(m, 'mp2'): # It is not a hairpin.
            m.p2 = getattr(*m.mp2)
            m.p3 = getattr(*m.mp3)
        m.p4 = getattr(*m.mp4)
        m.infer_points()
    return modules

def get_endpoints_arrowhead(modules):
    """
    
    Loops through modules and checks if any 5' ends connect to (touch) any 3' ends 
    and vice versa. If the answer is no, the 5' ends get a square, the 3' ends an arrow.
    Uses draw_arrowheads() from rendering.py.
    """
    
    #get endpoints
    fives = []
    threes = []
    
    for mod in modules:
    
        if len(mod.tentacle_data()[0]) == 3: #hairpin
            end5 = [mod.tentacle_data()[3][0]]
            end3 = [mod.tentacle_data()[4][2]]
        elif len(mod.tentacle_data()[0]) == 4: #4way
            end5 = [mod.tentacle_data()[3][0], mod.tentacle_data()[3][2]]
            end3 = [mod.tentacle_data()[4][3], mod.tentacle_data()[4][1]]
        else:
            print("bad module in get_endpoints")
            break
            
        fives.append(end5)
        threes.append(end3)
    
    #loop through and check closeness of opposite strand ends (except 0-0 and last-last)
    objects_inside = []
    
    curr_m = 0
    for m5 in fives:
        curr_end = 0
        for ends5 in m5:
            loose = True
            
            
            for m3 in threes:
                for ends3 in m3:
                    if almost_same(ends5, ends3):
                        loose = False
                
            if loose: #==True
                    
                i1, i2, i3, i4 = modules[curr_m].stem_data()[0:4]
                c_stem = modules[curr_m].stem_data()[4]
                    
                if curr_end == 0: #p1
                    if almost_same(fives[curr_m][curr_end], i1):
                        arrow5 = draw_arrowheads(i1, i2, c_stem, "5")
                    else:
                        c_tent = modules[curr_m].tentacle_data()[2][0] #p1 tent
                        c_tent = c_tent[0] #- first item
                            
                        arrow5 = draw_arrowheads(ends5, i1, c_tent, "5")
                        
                elif curr_end == 1: #p3
                    if almost_same(fives[curr_m][curr_end], i3):
                        arrow5 = draw_arrowheads(i3, i4, c_stem, "5")
                    else:
                        c_tent = modules[curr_m].tentacle_data()[2][2] #p3 tent
                        c_tent = c_tent[0] #- first item
                            
                        arrow5 = draw_arrowheads(ends5, i3, c_tent, "5")
                    
                objects_inside.append(arrow5)
                            
            curr_end += 1
        curr_m += 1
    
    #do it also with 3' ends
    curr_m = 0
    for m3 in threes:
        curr_end = 0
                
        for ends3 in m3:
            loose = True
            
            for m5 in fives:
                for ends5 in m5:
                    if almost_same(ends5, ends3):
                        loose = False        
                        
            if loose:

                i1, i2, i3, i4 = modules[curr_m].stem_data()[0:4]
                c_stem = modules[curr_m].stem_data()[4]
                    
                if curr_end == 0: #p4
                    if almost_same(threes[curr_m][curr_end], i4):
                        arrow3 = draw_arrowheads(i4, i3, c_stem, "3")
                    else:
                        c_tent = modules[curr_m].tentacle_data()[2][-1] #p4 tent
                        c_tent = c_tent[0] #- first item
                            
                        arrow3 = draw_arrowheads(ends3, i4, c_tent, "3")
                        
                elif curr_end == 1: #p2
                    if almost_same(threes[curr_m][curr_end], i2):
                        arrow3 = draw_arrowheads(i2, i1, c_stem, "3")
                    else:
                        c_tent = modules[curr_m].tentacle_data()[2][1] #p2 tent
                        c_tent = c_tent[0] #- first item
                            
                        arrow3 = draw_arrowheads(ends3, i2, c_tent, "3")
                    
                objects_inside.append(arrow3)
                
            
            curr_end += 1
        curr_m += 1
    
    return objects_inside

def push_apart_touching_ends(modules):
    """ If the main 5' and 3' ends touch, push them apart a visible amount."""
    
    first_p1 = modules[0].tentacle_data()[3][0]
    first_i1, first_i2 = modules[0].stem_data()[0:2]
    last_p4 = modules[-1].tentacle_data()[4][-1]
    last_i3, last_i4 = modules[-1].stem_data()[2:4]
    
    loose = True
    for mod in modules:
        if hasattr(mod, 'mp2'): # It is not a hairpin.
            p3 = mod.tentacle_data()[3][2]
            if almost_same(last_p4, p3):
                loose = False
    
    if loose: #if 3' on last module
        if almost_same(last_p4, first_p1): #if 5' and 3' touch
            #move these two points apart
            if almost_same(first_p1, first_i1): #no tentacle
                end_to_start = np.array([first_i1[0] - first_i2[0], first_i1[1] - first_i2[1]])
                unit_ets = end_to_start/np.linalg.norm(end_to_start)
                new_point = first_i1 - unit_ets*5
                
                modules[0].i1 = new_point #i1 rewritten
                modules[0].p1 = new_point #p1 rewritten

            else: #tentacle
                end_to_start = np.array([first_p1[0] - first_i1[0], first_p1[1] - first_i1[1]])
                unit_ets = end_to_start/np.linalg.norm(end_to_start)
                new_point = first_i1 - unit_ets*5
                
                modules[0].p1 = new_point #p1 rewritten
            
            #3'
            if almost_same(last_p4, last_i4): #no tentacle
                end_to_start = np.array([last_i4[0] - last_i3[0], last_i4[1] - last_i3[1]])
                unit_ets = end_to_start/np.linalg.norm(end_to_start)
                new_point = last_i4 - unit_ets*5
                
                modules[-1].i4 = new_point #i1 rewritten
                modules[-1].p4 = new_point #p1 rewritten

            else: #tentacle
                end_to_start = np.array([last_p4[0] - last_i4[0], last_p4[1] - last_i4[1]])
                unit_ets = end_to_start/np.linalg.norm(end_to_start)
                new_point = last_i4 - unit_ets*5
                
                modules[-1].p4 = new_point #p1 rewritten
        
    return modules


def draw_complex(cplx, pair_angles = None, loop_lengths = None, 
                 origin = (0, 0), scale = 10):
    """Returns the SVG image of a complex.

    The colors of domains have to be specified as property of the Domain
    objects: (DomainS.color = color). Pair_lengths and loop_angles are
    currently not supported, but may be usful extension for the future.

    Args:
        cplx (dsdobjects.Complex): A complex as given by dsd objects.
        pair_angles (list, optional): Provide the angles of paired domains in
                    order of their apparence in the specified complex rotation.
        loop_lengths (list, optional): Provide the end-to-end distance of loops
                    in order of their apparence in the specified complex rotation.

    Returns:
        svg: SVG image specification (drawsvg)
        pair_angles: the pair angles used to draw this svg image
        loop_lengths: the loop lengths used to draw this svg image

    """
    pa, ll = get_default_plot_params(cplx)
    assert pair_angles is None or len(pair_angles) == len(pa)
    assert loop_lengths is None or len(loop_lengths) == len(ll)

    if pair_angles is None:
        pair_angles = pa
    if loop_lengths is None:
        loop_lengths = ll

    modules = get_final_modules(cplx, pair_angles, loop_lengths, 
                                origin = origin, scale = scale)

    modules = push_apart_touching_ends(modules)

    objects = []
    
    for m in modules:
        objects.extend(draw_stem(*m.stem_data()))
        [objects.extend(x) for x in draw_tentacles(*m.tentacle_data())]
    
    arrows = get_endpoints_arrowhead(modules)
    objects.extend(arrows)
                
    return objects, pair_angles, loop_lengths