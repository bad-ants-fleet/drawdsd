
import dsdobjects
from dsdobjects import DomainS
from dsdobjects.objectio import set_io_objects, clear_io_objects, read_pil, read_pil_line
from itertools import combinations

import rendering
import drawmodules as dms

def agl(a):
    return a % 360

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
                m = dms.dsd_hairpin_module(stem, t1, t4, th, scale = scale)
                # Provide the loop_length data
                if k1: m.k1 = scale * k1
                if k4: m.k4 = scale * k4
                # We pretend it has four ends to save some code later.
                meta = [[*sg1[0], [si, di], *sg1[1]], 
                        [[sj, dj], *sg2[1]]]
            else:
                t2 = [stable[i][j] for (i, j) in sg1[1]]
                t3 = [stable[i][j] for (i, j) in sg2[0]]
                m = dms.dsd_module(stem, t1, t2, t3, t4, scale = scale)
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

def draw_complex(cplx, pair_angles = None, loop_lengths = None):
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

    modules = get_final_modules(cplx, pair_angles, loop_lengths, origin = (0, 0), scale = 10)

    svg = rendering.get_canvas()
    for m in modules:
        svg.extend(rendering.draw_stem(*m.stem_data()))
        [svg.extend(x) for x in rendering.draw_tentacles(*m.tentacle_data())]

    return svg, pair_angles, loop_lengths

def main():
    dcolors = {'r': 'red',
           'r*': 'red',
           'a': 'blue',
           'a*': 'blue',
           'b': 'blue',
           'b*': 'blue',
           'c': 'blue',
           'c*': 'blue',
           'g': 'green',
           'g*': 'green',
           'x': 'green',
           'x*': 'green',
           'k': 'green',
           'k*': 'green',
           'l': 'black',
           'l*': 'black',
           'y': 'yellow',
           'y*': 'yellow'}

    domains = '''
    length a = 7
    length b = 5
    length c = 15
    length x = 5
    length y = 15
    length k = 5
    length r = 5
    length g = 5
    length l = 5
    '''
    set_io_objects()
    _ = read_pil(domains)

    # Test some examples 
    #mycplx = read_pil_line('A = k a( b y y c( y + ) ) k')
    mycplx = read_pil_line('A = r b( g r b( l ) y r b( g + l ) y r b( g + l ) y l ) y')
    #mycplx = read_pil_line('A = r b( b( l ) b( g + l ) l ) y')

    for dom in mycplx.domains:
        dom.color = dcolors[dom.name]

    # provide module angles in order
    #'A = x <@0> a( <|2|> b y y <@0> c( y + ) ) k'
    # change default tentacle lengths

    # This assumes you provide a complex e.g. from dsdobjects, and also have
    # colors assigned to domains.
    pa, ll = get_default_plot_params(mycplx)

    # mess with defaults:
    pa[3] = 180
    ll[0][1] = 20
    ll[1][1] = 20

    svg, pa, ll = draw_complex(mycplx, pair_angles = pa, loop_lengths = ll)

    svg.savePng('example.png')
    #print(d.asSvg())

if __name__ == '__main__':
    main()

