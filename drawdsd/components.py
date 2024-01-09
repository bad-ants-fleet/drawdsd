#
# Provides "components" (or "drawingmodules") to calculate plot coordinates.
#
import logging
log = logging.getLogger(__name__)

import numpy as np

class DrawingModuleError(Exception):
    pass

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

class fourway_module:
    def __init__(self, pair, t1, t2, t3, t4, p15, p23, p35, p43, scale = 10):
        self.pair = pair
        assert len(pair) == 2
        assert len(pair[0]) == len(pair[1])
        self.plength = scale * (len(pair[0])+1)
        self.pnameT = pair[0].name
        self.pnameB = pair[1].name
        assert pair[0].color == pair[1].color
        self.pcolor = pair[0].color

        self.p15 = p15
        self.p23 = p23
        self.p35 = p35
        self.p43 = p43

        self.n1 = [dom.name for dom in t1] 
        self.l1 = [scale * (dom.length+1) for dom in t1] 
        self.k1 = sum(self.l1)
        self.c1 = [dom.color for dom in t1] 
        self.m1 = [None for dom in t1] 
        if t1 and p15: 
            self.m1[0] = 'p5' 
            self.p15 = False

        self.n2 = [dom.name for dom in t2] 
        self.l2 = [scale * (dom.length+1) for dom in t2] 
        self.k2 = sum(self.l2)
        self.c2 = [dom.color for dom in t2] 
        self.m2 = [False for dom in t2] 
        if t2 and p23: 
            self.m2[0] = 'p3' 
            self.p23 = False

        self.n3 = [dom.name for dom in t3] 
        self.l3 = [scale * (dom.length+1) for dom in t3] 
        self.k3 = sum(self.l3)
        self.c3 = [dom.color for dom in t3] 
        self.m3 = [None for dom in t3] 
        if t3 and p35:
            self.m3[0] = 'p5' 
            self.p35 = False

        self.n4 = [dom.name for dom in t4] 
        self.l4 = [scale * (dom.length+1) for dom in t4] 
        self.k4 = sum(self.l4)
        self.c4 = [dom.color for dom in t4] 
        self.m4 = [False for dom in t4] 
        if t4 and p43: 
            self.m4[0] = 'p3' 
            self.p43 = False
        

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

        TODO: Need implementation for a complex where p1 is None!
            'A = x( r + r( + ) x*( + ) t2*( + ) )'
        """
        if self.p1 is None:
            raise NotImplementedError('No starting point given!')
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

    def stem_data(self):
        return (self.i1, self.i2, self.i3, self.i4, 
                self.p15, self.p23, self.p35, self.p43, 
                self.pcolor, self.angle, self.plength, 
                self.pnameT, self.pnameB)

    def tentacle_data(self):
        # paired domains
        dns = [self.n1, self.n2, self.n3, self.n4]
        dls = [self.l1, self.l2, self.l3, self.l4]
        dcs = [self.c1, self.c2, self.c3, self.c4]
        dzs = [self.p1, self.i2, self.p3, self.i4]
        dos = [self.i1, self.p2, self.i3, self.p4]
        dms = [self.m1, self.m2, self.m3, self.m4]
        return dns, dls, dcs, dzs, dos, dms

class hairpin_module:
    def __init__(self, pair, t1, t4, th, p15, p43, scale = 10):
        self.pair = pair
        assert len(pair) == 2
        assert len(pair[0]) == len(pair[1])
        self.plength = scale * (len(pair[0])+1)
        self.pnameT = pair[0].name
        self.pnameB = pair[1].name
        assert pair[0].color == pair[1].color
        self.pcolor = pair[0].color

        self.p15 = p15
        self.p43 = p43


        self.n1 = [dom.name for dom in t1] 
        self.l1 = [scale * (dom.length+1) for dom in t1] 
        self.k1 = sum(self.l1)
        self.c1 = [dom.color for dom in t1] 
        self.m1 = [False for dom in t1] 
        if t1 and p15: 
            self.m1[0] = 'p5' 
            self.p15 = False
        self.n4 = [dom.name for dom in t4] 
        self.l4 = [scale * (dom.length+1) for dom in t4] 
        self.k4 = sum(self.l4)
        self.c4 = [dom.color for dom in t4] 
        self.m4 = [False for dom in t4] 
        if t4 and p43: 
            self.m4[0] = 'p3' 
            self.p43 = False
        self.nh = [dom.name for dom in th] 
        self.lh = [scale * (dom.length+1) for dom in th] 
        self.ch = [dom.color for dom in th] 
        self.mh = [False for dom in th] 

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

    def stem_data(self):
        return (self.i1, self.i2, self.i3, self.i4, 
                self.p15, False, False, self.p43, 
                self.pcolor, self.angle, self.plength, 
                self.pnameT, self.pnameB)

    def tentacle_data(self):
        # paired domains
        dns = [self.n1, self.nh, self.n4]
        dls = [self.l1, self.lh, self.l4]
        dcs = [self.c1, self.ch, self.c4]
        dzs = [self.p1, self.i2, self.i4]
        dos = [self.i1, self.i3, self.p4]
        dms = [self.m1, self.mh, self.m4]
        return dns, dls, dcs, dzs, dos, dms

