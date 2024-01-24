#
# High-level and mid-level library interface.
#
from itertools import combinations

from .components import fourway_module, hairpin_module
from .rendering import get_drawing, get_rgb_palette, draw_stem, draw_tentacles

def agl(a):
    return a % 360

def get_default_plot_params(stable, ptable):
    """ Produce the defaults for:
        - pair_angles
        - loop_lengths
    """
    pair_angles = []
    loop_lengths = []
    loop_angles = []
    for si, strand in enumerate(ptable):
        nlen, nang, llen = [], [], 0
        for di, pair in enumerate(strand):
            if pair is None: # unpaired
                llen += (stable[si][di].length)
            else: # pair
                if pair > (si, di):
                    pair_angles.append(0)
                nlen.append(llen)
                nang.append(None)
                llen = 0
        nlen.append(llen)
        nang.append(None)
        loop_lengths.append(nlen)
        loop_angles.append(nang)
    return pair_angles, loop_lengths, loop_angles

def get_segments(ptable, loop_lengths, loop_angles):
    """Helper function to decompose a pair table into segments.

    Segments[(pair, pair)] = ([segment1, segment2], (ll1, ll2), (5' end, 3' end)),

    """
    segments = {}
    identity, seg, pos = None, [[], []], 0
    ll1, ll2, lli = 0, 0, 0 # loop_length stuff
    p5, p3 = True, False
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pos:
                ll2 = loop_lengths[si][lli]
                la2 = loop_angles[si][lli]
            else:
                ll1 = loop_lengths[si][lli]
                la1 = loop_angles[si][lli]
            if pair is not None:
                pos += 1
                lli += 1
                if identity is not None:
                    segments[identity] = seg, (ll1, ll2), (la1, la2), (p5, p3)
                    identity, seg, pos = None, [[], []], 1
                    ll1, ll2 = 0, 0
                    p5, p3 = False, False
                identity = (si, di), pair
            else:
                seg[pos].append([si, di])
        # Strand break means we start a new segment!
        p3 = True
        segments[identity] = seg, (ll1, ll2), (la1, la2), (p5, p3)
        identity, seg, pos = None, [[], []], 0
        ll1, ll2, lli = 0, 0, 0
        p5, p3 = True, False
    return segments

def get_raw_modules(stable, ptable, segments, pair_angles):
    """Initialize the 4-way and the hairpin module.
    """
    modules, metadata = [], []
    # We start at 5' end, so the modules must be in order!
    pa = 0
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pair is None or (si, di) > pair:
                continue
            (sj, dj) = pair
            # Get both segments by their identifier.
            sg1, (k1, k2), (fa1, fa2), (p15, p23) = segments[(si, di), (sj, dj)]
            sg2, (k3, k4), (fa3, fa4), (p35, p43) = segments[(sj, dj), (si, di)]
            # The stem is formed via the single paired domain.
            stem = [stable[si][di], stable[sj][dj]] 
            # The tentacles listed in segments.
            t1 = [stable[i][j] for (i, j) in sg1[0]]
            t4 = [stable[i][j] for (i, j) in sg2[1]]
            # Check if it is a hairpin!
            if (sg1[1] and sg2[0] == [] and sg1[1][-1] == [sj, dj-1]) or \
                    (sg1[1] == [] and sg2[0] == [] and (si, di+1) == (sj, dj) and not p23 and not p35):
                th = [stable[i][j] for (i, j) in sg1[1]]
                m = hairpin_module(stem, t1, t4, th, p15, p43)
                # Provide the loop_length data
                if k1: m.k1 = k1
                if k4: m.k4 = k4
                m.fa1 = fa1 if fa1 is None else -fa1 % 360
                m.fa4 = fa4 if fa4 is None else -fa4 % 360
                # We pretend it has four ends to save some code later.
                meta = [[*sg1[0], [si, di], *sg1[1]], 
                        [         [sj, dj], *sg2[1]]]
            else:
                t2 = [stable[i][j] for (i, j) in sg1[1]]
                t3 = [stable[i][j] for (i, j) in sg2[0]]
                m = fourway_module(stem, t1, t2, t3, t4, p15, p23, p35, p43)
                # Provide the loop_length data
                if k1: m.k1 = k1
                if k2: m.k2 = k2
                if k3: m.k3 = k3
                if k4: m.k4 = k4
                m.fa1 = fa1 if fa1 is None else -fa1 % 360
                m.fa2 = fa2 if fa2 is None else -fa2 % 360
                m.fa3 = fa3 if fa3 is None else -fa3 % 360
                m.fa4 = fa4 if fa4 is None else -fa4 % 360
                meta = [[*sg1[0], [si, di], *sg1[1]], 
                        [*sg2[0], [sj, dj], *sg2[1]]]
            # Provide the pair_angle data
            m.angle = -pair_angles[pa]
            pa += 1 
            modules.append(m)
            metadata.append(meta)
    return modules, metadata

def get_final_modules(stable, ptable, pair_angles, loop_lengths, loop_angles, origin = (0, 0)):

    # Get the helper datastructure segments.
    segments = get_segments(ptable, loop_lengths, loop_angles)

    # Returns all modules but has no knowledge on their neighborhood.
    modules, metadata = get_raw_modules(stable, ptable, segments, pair_angles)

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

    # Set the angles for all loop regions connecting stems.
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
    done = []
    maxiter = len(modules) * len(modules)
    it = 0
    while len(modules):
        m = modules.pop(0)
        m.angle = m.angle # Re-setting the angle updates all small angles!
        m.p1 = getattr(*m.mp1)
        if hasattr(m, 'mp2'): # It is not a hairpin.
            m.p2 = getattr(*m.mp2)
            m.p3 = getattr(*m.mp3)
        m.p4 = getattr(*m.mp4)
        if m.infer_points_p1():
            done.append(m)
        elif m.infer_points_p4():
            done.append(m)
        else:
            modules.append(m)
        if it > maxiter:
            raise SystemExit('Disconnected Complexes')
        it += 1 
    modules = done
    return modules

def get_svg_components(stable, ptable, pair_angles = None, loop_lengths = None, loop_angles = None,
                 rotate = 0, spacing = 0, origin = (0, 0)):
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
    pa, ll, la = get_default_plot_params(stable, ptable)
    assert pair_angles is None or len(pair_angles) == len(pa)
    assert loop_lengths is None or len(loop_lengths) == len(ll)
    assert loop_angles is None or len(loop_angles) == len(la)

    if pair_angles is None:
        pair_angles = pa
    if loop_lengths is None:
        loop_lengths = ll
    if loop_angles is None:
        loop_angles = la

    # Provide spacing
    loop_lengths = [[(i+spacing) if (i==0) else i for i in s] for s in loop_lengths]
    # Rotate complex
    pair_angles = [a+rotate for a in pair_angles]

    modules = get_final_modules(stable, ptable, pair_angles, loop_lengths, 
                                loop_angles,
                                origin = origin)
    objects = []
    for m in modules:
        objects.extend(draw_stem(*m.stem_data()))
        [objects.extend(x) for x in draw_tentacles(*m.tentacle_data())]

    return objects, pair_angles, loop_lengths, loop_angles

