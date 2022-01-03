
import dsdobjects
from itertools import combinations
from dsdobjects import DomainS
from dsdobjects.objectio import set_io_objects, clear_io_objects, read_pil, read_pil_line

import drawSvg as draw
import numpy as np

class DrawingError(Exception):
    pass

dfont = 'bold'
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


def agl(a):
    return a % 360

def half_angle(a, b, fwd = False, rev = False):
    # from a -> to b
    assert b == (b % 360)
    assert a == (a % 360)
    if a != 0 and b == 0 and fwd:
        b = 360 
    if b != 0 and a == 0 and rev:
        a = 360 
    return (b-a)/2

def get_coords(i1, ang, c = 1):
    assert ang == (ang % 360)
    x0, y0 = i1
    if ang <= 90:
        rad = np.radians(ang)
        return x0+(np.cos(rad)*c), y0+(np.sin(rad)*c)
    elif 90 < ang <= 180:
        rad = np.radians(180 - ang)
        return x0-(np.cos(rad)*c), y0+(np.sin(rad)*c)
    elif 180 < ang <= 270:
        rad = np.radians(ang - 180)
        return x0-(np.cos(rad)*c), y0-(np.sin(rad)*c)
    else:
        rad = np.radians(360 - ang)
        return x0+(np.cos(rad)*c), y0-(np.sin(rad)*c)

def get_a1(pa, a): # red, from pa to a
    """ Get angle from i1 to p1.
    if pa=180 and a=0, then ang=270!
    (previous is below current!)
    """
    #print('a1', pa, a)
    if pa is None:
        return agl(a-30+180) # from i1
    rel = agl(pa-a)
    if 0 < rel < 180:
        return agl(a + rel/2 + 180)
    if 180 <= rel < 360:
        return agl(a + rel/2)
    assert np.isclose(rel, 0)
    return agl(a + 180)

def get_a2(a, na): # green, from a to na
    #print('a2', a, na)
    if na is None:
        return agl(a+30) # from i2
    rel = agl(na-a)
    if 0 < rel <= 180:
        return agl(a + rel/2)
    if 180 < rel < 360:
        return agl(a + rel/2 + 180)
    assert np.isclose(rel, 0)
    return a

def get_a3(pa, a): # black, from a to na
    if pa is None:
        return agl(a-30+180) # from i3

    rel = agl(pa-a)
    #print('a3', pa, a, a-pa, rel)
    if 0 < rel < 180:
        return agl(a + rel/2 + 180)
    if 180 <= rel < 360:
        return agl(a + rel/2)
    assert np.isclose(rel, 0)
    return agl(a + 180)

def get_a4(a, na): # yellow, from a to na
    #print('a4', a, na)
    if na is None:
        return agl(a+30) # from i4
    """
    180 -> 0 == 0 -> 180 

    """
    rel = agl(na-a)
    if 0 < rel <= 180:
        return agl(a + rel/2)
    if 180 < rel < 360:
        return agl(a + rel/2 + 180)
    assert np.isclose(rel, 0)
    return agl(a)

def get_ltable(cplx):
    stable = list(cplx.strand_table)
    ptable = list(cplx.pair_table)
    ltable = []
    for si, strand in enumerate(ptable):
        slen = []
        l = 0
        for di, pair in enumerate(strand):
            if pair is None:
                l += stable[si][di].length
            else:
                slen.append(l)
                l = 0
        slen.append(l)
        ltable.append(slen)
    return ltable
 

def decompose(cplx, ltable):
    # print(cplx)
    # print(cplx.kernel_string)
    # print(list(cplx.sequence))
    # print(list(cplx.structure))
    # print(list(cplx.strand_table))
    # print(list(cplx.pair_table))
    stable = list(cplx.strand_table)
    ptable = list(cplx.pair_table)

    segments = {}
    identity, seg, pos = None, [[], []], 0
    l1, l2, li = 0, 0, 0
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pair is not None:
                pos += 1
                li += 1
                if identity is not None:
                    segments[identity] = seg, (l1, l2)
                    identity, seg, pos = None, [[], []], 1
                    l1, l2 = 0, 0
                identity = (si, di), pair
            else:
                if pos:
                    l2 = ltable[si][li]
                else:
                    l1 = ltable[si][li]
                seg[pos].append([si, di])
        # Strand break means we start a new segment!
        segments[identity] = seg, (l1, l2)
        identity, seg, pos = None, [[], []], 0
        l1, l2, li = 0, 0, 0
    return segments

def get_modules(cplx, segments, angles, scale):
    modules = []
    meta = []
    stable = list(cplx.strand_table)
    ptable = list(cplx.pair_table)

    # We start at ptable 5' end!
    for si, strand in enumerate(ptable):
        for di, pair in enumerate(strand):
            if pair is None:
                continue
            (sj, dj) = pair
            if (si, di) > (sj, dj):
                # No double counting.
                continue
            # Get both segments by their identifier.
            sg1, (k1, k2) = segments[(si, di), (sj, dj)]
            sg2, (k3, k4) = segments[(sj, dj), (si, di)]

            # The stem is formed via the single paired domain.
            ST = [stable[si][di], stable[sj][dj]] 
            D1 = [stable[i][j] for (i, j) in sg1[0]]
            D4 = [stable[i][j] for (i, j) in sg2[1]]
            # Check if it is a hairpin!
            if sg1[1] and sg2[0] == [] and sg1[1][-1] == [sj, dj-1]:
                H = [stable[i][j] for (i, j) in sg1[1]]
                m = dsd_hairpin_module(ST, D1, D4, H, scale = scale)
                if k1: m.k1 = scale * k1
                if k4: m.k4 = scale * k4
                mt = [[*sg1[0], [si, di], *sg1[1]], [[sj, dj], *sg2[1]]]
            else:
                D2 = [stable[i][j] for (i, j) in sg1[1]]
                D3 = [stable[i][j] for (i, j) in sg2[0]]
                m = dsd_module(ST, D1, D2, D3, D4, scale = scale)
                if k1: m.k1 = scale * k1
                if k2: m.k2 = scale * k2
                if k3: m.k3 = scale * k3
                if k4: m.k4 = scale * k4
                mt = [[*sg1[0], [si, di], *sg1[1]], 
                      [*sg2[0], [sj, dj], *sg2[1]]]
            m.angle = angles.pop(0)
            modules.append(m)
            meta.append(mt)

    # Set all the small angles
    for (m, mt), (n, nt) in combinations(zip(modules, meta), 2):
        md1, md2 = mt[0][0], mt[0][-1] # first and last
        md3, md4 = mt[1][0], mt[1][-1]
        nd1, nd2 = nt[0][0], nt[0][-1]
        nd3, nd4 = nt[1][0], nt[1][-1]

        # The strand continues from m(2) to n(1)
        if (md2[0], md2[1]) == (nd1[0], nd1[1]-1):
            m._nma2 = n.angle
            n._pma1 = m.angle
            assert n.mp1[0] == n
            assert m.mp2[0] == m
            n.mp1 = (m, 'p2')
            m.mp2 = (n, 'p1')

        # The strand continues from n(2) to m(1)
        if (nd2[0], nd2[1]) == (md1[0], md1[1]-1):
            n._nma2 = m.angle
            m._pma1 = n.angle
            assert m.mp1[0] == m
            assert n.mp2[0] == n
            m.mp1 = (n, 'p2')
            n.mp2 = (m, 'p1')

        # The strand continues from m(4) to n(1)
        if (md4[0], md4[1]) == (nd1[0], nd1[1]-1):
            m._nma4 = n.angle
            n._pma1 = agl(m.angle + 180)
            assert n.mp1[0] == n
            assert m.mp4[0] == m
            n.mp1 = (m, 'p4')
            m.mp4 = (n, 'p1')

        # The strand continues from n(4) to m(1)
        if (nd4[0], nd4[1]) == (md1[0], md1[1]-1):
            n._nma4 = m.angle
            n._pma1 = agl(n.angle + 180)
            assert m.mp1[0] == m
            assert n.mp4[0] == n
            m.mp1 = (n, 'p4')
            n.mp4 = (m, 'p1')

        # The strand continues from m(4) to n(3)
        if (md4[0], md4[1]) == (nd3[0], nd3[1]-1):
            m._nma4 = agl(n.angle + 180)
            n._pma3 = agl(m.angle + 180)
            assert n.mp3[0] == n
            assert m.mp4[0] == m
            n.mp3 = (m, 'p4')
            m.mp4 = (n, 'p3')

        # The strand continues from n(4) to m(3)
        if (nd4[0], nd4[1]) == (md3[0], md3[1]-1):
            n._nma4 = agl(m.angle + 180)
            m._pma3 = agl(n.angle + 180)
            assert m.mp3[0] == m
            assert n.mp4[0] == n
            m.mp3 = (n, 'p4')
            n.mp4 = (m, 'p3')

        # The strand continues from n(2) to m(3)
        if (nd2[0], nd2[1]) == (md3[0], md3[1]-1):
            raise NotImplementedError('needed?')
        # The strand continues from m(2) to n(3)
        if (md2[0], md2[1]) == (nd3[0], nd3[1]-1):
            raise NotImplementedError('needed?')


    m = modules[0]
    m.p1 = (0, 0)
    for m in modules:
        m.angle = m.angle
        m.p1 = getattr(*m.mp1)
        if hasattr(m, 'mp2'):
            m.p2 = getattr(*m.mp2)
            m.p3 = getattr(*m.mp3)
        m.p4 = getattr(*m.mp4)
        m.infer_points()
    return modules

def draw_center(m): # Takes a module!
    bground = draw.Path(stroke_width=0, fill='black', fill_opacity=.3)
    dompath = draw.Path(stroke_width=5, stroke=m.pcolor, fill_opacity=0)

    (x1, y1) = m.i1
    (x2, y2) = m.i2
    (x3, y3) = m.i3
    (x4, y4) = m.i4

    bground.m(x1, y1).L(x2, y2).L(x3, y3).L(x4, y4)
    dompath.m(x1, y1).L(x2, y2).M(x3, y3).L(x4, y4)

    # TODO ... dirty hack
    if m.angle < 45:
        (cT, cB) = (10, 20)
    elif 45 <= m.angle <= 135:
        (cT, cB) = (10, 15)
    elif 135 < m.angle < 225:
        (cT, cB) = (20, 10)
    else:
        (cT, cB) = (10, 15)

    tx, ty = (x1+x2)/2, (y1+y2)/2 # point in the middle on straight line
    tx, ty = tx-(y2-y1)/m.plength*cT, ty+(x2-x1)/m.plength*cT
    tagT = draw.Text(m.pnameT, 14, x=tx, y=ty, font_weight = dfont, text_anchor='middle', valign='center')

    tx, ty = (x3+x4)/2, (y3+y4)/2 # point in the middle on straight line
    tx, ty = tx-(y4-y3)/m.plength*cB, ty+(x4-x3)/m.plength*cB
    tagB = draw.Text(m.pnameB, 14, x=tx, y=ty, font_weight = dfont, text_anchor='middle', valign='center')

    return bground, dompath, tagT, tagB

def draw_tentacles(dns, dls, dcs, dzs, dos):
    # paired domains
    for ns, ls, cs, (x0, y0), (x1, y1) in zip(dns, dls, dcs, dzs, dos):
        num, length = len(ns), sum(ls)
        linelen = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        if not np.isclose(length, linelen) and linelen > length:
            ns.append('')
            ls.append(linelen-length)
            cs.append('black')
            length = sum(ls)
            num += 1 
        if np.isclose(length, linelen):
            xs, ys = x0, y0
            for i in range(num):
                xn, yn = xs+(x1-x0)/length*ls[i], ys+(y1-y0)/length*ls[i]

                if ns[i] == '':
                    da = ' '.join(["5" for _ in range(int(round(ls[i]+1)))])
                    p = draw.Path(stroke_width=2, stroke=cs[i], fill_opacity=0, 
                        stroke_dasharray=da)
                else:
                    p = draw.Path(stroke_width=5, stroke=cs[i], fill_opacity=0)
                p.m(xs, ys).L(xn, yn)
                cT = 20 if xs > xn else 10
                tx, ty = (xs+xn)/2, (ys+yn)/2 # point in the middle on straight line
                tx, ty = tx-(yn-ys)/ls[i]*cT, ty+(xn-xs)/ls[i]*cT
                t = draw.Text(ns[i], 14, x=tx, y=ty, font_weight = dfont, 
                              text_anchor='middle', valign='center')
                yield p, t
                xs, ys = xn, yn
        elif length > linelen:
            missing = length - linelen
            # The points of the rectangle: M(x0, y0).L(rx2, ry2).L(rx3, ry3).L(x1, y1)
            (rx2, ry2) = x0-(y1-y0)/linelen*missing/2, y0+(x1-x0)/linelen*missing/2
            (rx3, ry3) = x1-(y1-y0)/linelen*missing/2, y1+(x1-x0)/linelen*missing/2

            llen = linelen
            slen = np.sqrt((rx2-x0)**2 + (ry2-y0)**2)
            assert np.isclose(length, llen+2*slen)

            xs, ys = x0, y0
            cseg, ci = [(rx2, ry2), (rx3, ry3), (x1, y1)], 0
            for i in range(num):
                l = ls[i]
                tl = ls[i]/2# + 5
                p = draw.Path(stroke_width=5, stroke=cs[i], fill_opacity=0)
                p.M(xs, ys)
                while True:
                    xn, yn = cseg[ci]
                    dlen = np.sqrt((xn-xs)**2 + (yn-ys)**2)
                    if not np.isclose(dlen, tl) and dlen < tl:
                        tl -= dlen
                    else:
                        tx, ty = xs+(xn-xs)/dlen*tl, ys+(yn-ys)/dlen*tl

                        cT = 20 if xs > xn else 10
                        tx, ty = tx-(yn-ys)/dlen*cT, ty+(xn-xs)/dlen*cT
                        t = draw.Text(ns[i], 14, x=tx, y=ty, font_weight = dfont,
                                      text_anchor='middle', valign='center')
                        tl = 2 * l
                    if not np.isclose(dlen, l) and dlen < l:
                        p.L(xn, yn)#.L(rx3, ry3).L(x1, y1)
                        l -= dlen
                        ci += 1 
                        xs, ys = xn, yn
                    else:
                        xn, yn = xs+(xn-xs)/dlen*l, ys+(yn-ys)/dlen*l
                        p.L(xn, yn)#.L(rx3, ry3).L(x1, y1)
                        xs, ys = xn, yn
                        break
                yield p, t
        else:
            raise NotImplementedError(linelen, length)

class dsd_module:
    def __init__(self, pair, t1, t2, t3, t4, scale = 10):

        self.pair = pair
        assert len(pair) == 2
        assert len(pair[0]) == len(pair[1])
        self.plength = scale * len(pair[0])
        self.pnameT = pair[0].name
        self.pnameB = pair[1].name
        assert dcolors[pair[0].name] == dcolors[pair[1].name]
        self.pcolor = dcolors[self.pnameT]

        self.n1 = [dom.name for dom in t1] 
        self.l1 = [scale * dom.length for dom in t1] 
        self.k1 = sum(self.l1)
        self.c1 = [dcolors[dom.name] for dom in t1] 
        self.n2 = [dom.name for dom in t2] 
        self.l2 = [scale * dom.length for dom in t2] 
        self.k2 = sum(self.l2)
        self.c2 = [dcolors[dom.name] for dom in t2] 
        self.n3 = [dom.name for dom in t3] 
        self.l3 = [scale * dom.length for dom in t3] 
        self.k3 = sum(self.l3)
        self.c3 = [dcolors[dom.name] for dom in t3] 
        self.n4 = [dom.name for dom in t4] 
        self.l4 = [scale * dom.length for dom in t4] 
        self.k4 = sum(self.l4)
        self.c4 = [dcolors[dom.name] for dom in t4] 

        # get/set angles
        self._pma1 = None # prev module angle 1
        self._nma2 = None # next module angle 2
        self._pma3 = None
        self._nma4 = None
        self.angle = 0 # Initializes self.a1, etc.

        # get/set exterior points 
        self.p1 = None
        self.p2 = None
        self.p3 = None
        self.p4 = None

        # initialize point cross-referencing 
        # TODO: is that a good way to do it?
        self.mp1 = (self, 'p1')
        self.mp2 = (self, 'p2') 
        self.mp3 = (self, 'p3') 
        self.mp4 = (self, 'p4') 
 
        # initialize internal points 
        self.i1 = None
        self.i2 = None
        self.i3 = None
        self.i4 = None

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle):
        """Sets the angle, and adjusts all small angles accordingly."""
        self._angle = agl(angle)
        self.a1 = get_a1(self._pma1, self.angle)
        self.a2 = get_a2(self.angle, self._nma2)
        self.a3 = get_a3(self._pma3, agl(self.angle+180))
        self.a4 = get_a4(agl(self.angle+180), self._nma4)

    def infer_points(self):
        """Infer remaining points given p1 and the angles.

        TODO: Is there a situation where we need a different
        starting point than p1?
        """
        if self.p1 is None:
            raise DrawingError('No starting point given!')
        if self.i1 is None: # Reverse direction (+180 angle)
            self.i1 = get_coords(self.p1, agl(self.a1+180), self.k1)
        if self.i2 is None:
            self.i2 = get_coords(self.i1, self.angle, self.plength)
        if self.i3 is None:
            self.i3 = get_coords(self.i2, agl(self.angle+270), 25)
        if self.i4 is None:
            self.i4 = get_coords(self.i1, agl(self.angle+270), 25)
        if self.p2 is None:
            self.p2 = get_coords(self.i2, self.a2, self.k2)
        if self.p3 is None:
            self.p3 = get_coords(self.i3, self.a3, self.k3)
        if self.p4 is None:
            self.p4 = get_coords(self.i4, self.a4, self.k4)

    def draw_center(self):
        return draw_center(self)

    def draw_tentacles(self):
        # paired domains
        dns = [self.n1, self.n2, self.n3, self.n4]
        dls = [self.l1, self.l2, self.l3, self.l4]
        dcs = [self.c1, self.c2, self.c3, self.c4]
        dzs = [self.p1, self.i2, self.p3, self.i4]
        dos = [self.i1, self.p2, self.i3, self.p4]
        for x in draw_tentacles(dns, dls, dcs, dzs, dos):
            yield x

class dsd_hairpin_module:
    def __init__(self, pair, t1, t4, th, scale = 10):
        self.pair = pair
        assert len(pair) == 2
        assert len(pair[0]) == len(pair[1])
        self.plength = scale * len(pair[0])
        self.pnameT = pair[0].name
        self.pnameB = pair[1].name
        assert dcolors[pair[0].name] == dcolors[pair[1].name]
        self.pcolor = dcolors[self.pnameT]

        self.n1 = [dom.name for dom in t1] 
        self.l1 = [scale * dom.length for dom in t1] 
        self.k1 = sum(self.l1)
        self.c1 = [dcolors[dom.name] for dom in t1] 
        self.n4 = [dom.name for dom in t4] 
        self.l4 = [scale * dom.length for dom in t4] 
        self.k4 = sum(self.l4)
        self.c4 = [dcolors[dom.name] for dom in t4] 
        self.nh = [dom.name for dom in th] 
        self.lh = [scale * dom.length for dom in th] 
        self.ch = [dcolors[dom.name] for dom in th] 

        # get/set Angle/Coordinates
        self._pma1 = None
        self._nma4 = None
        self.angle = 0

        self.p1 = None
        self.p4 = None
        self.mp1 = (self, 'p1')
        self.mp4 = (self, 'p4')

        self.i1 = None
        self.i2 = None
        self.i3 = None
        self.i4 = None

    @property
    def angle(self):
        return self._angle

    @angle.setter
    def angle(self, angle):
        self._angle = agl(angle)
        self.a1 = get_a1(self._pma1, self.angle)
        self.a4 = get_a4(agl(self.angle+180), self._nma4)

    def infer_points(self):
        # autocomplete all coordinates.
        if self.i1 is None: # Reverse dirction (+180 angle)
            self.i1 = get_coords(self.p1, agl(self.a1 + 180), self.k1)
        if self.i2 is None:
            self.i2 = get_coords(self.i1, self.angle, self.plength)
        if self.i3 is None:
            self.i3 = get_coords(self.i2, agl(self.angle+270), 25)
        if self.i4 is None:
            self.i4 = get_coords(self.i1, agl(self.angle+270), 25)
        if self.p4 is None:
            self.p4 = get_coords(self.i4, self.a4, self.k4)

    def draw_center(self):
        return draw_center(self)

    def draw_tentacles(self):
        # paired domains
        dns = [self.n1, self.nh, self.n4]
        dls = [self.l1, self.lh, self.l4]
        dcs = [self.c1, self.ch, self.c4]
        dzs = [self.p1, self.i2, self.i4]
        dos = [self.i1, self.i3, self.p4]
        for x in draw_tentacles(dns, dls, dcs, dzs, dos):
            yield x

def main():
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

    # provide module angles in order
    #'A = x <@0> a( <|2|> b y y <@0> c( y + ) ) k'
    angles = [0, 0, 0, 180]

    # change default tentacle lengths
    ltable = get_ltable(mycplx)
    #print(ltable)
    ltable[0][1] = 5
    ltable[1][1] = 5

    segments = decompose(mycplx, ltable)
    ms = get_modules(mycplx, segments, angles, scale = 10)

    # Initialize the SVG image & background
    d = draw.Drawing(800, 800, origin = (-250, -550), displayInline = False)
    d.append(draw.Rectangle(-100, -400, 500, 500, fill='white'))

    for m in ms:
        d.extend(m.draw_center())
        [d.extend(x) for x in m.draw_tentacles()]

    d.savePng('example.png')
    #print(d.asSvg())

if __name__ == '__main__':
    main()

