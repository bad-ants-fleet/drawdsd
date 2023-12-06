#!/usr/bin/env python

import unittest
import numpy as np

from drawdsd.components import half_angle, get_coords, get_a1, get_a2, get_a3, get_a4

SKIP = False

class TestUtils(unittest.TestCase):
    def test_half_angle(self):
        assert   0 == half_angle(0, 0)
        assert   0 == half_angle(0, 0, fwd = True)
        #assert   5 == half_angle(350, 0)
        assert   5 == half_angle(350, 0, fwd = True)

        assert  45 == half_angle(  0,  90)
        assert  90 == half_angle(  0, 180)
        assert 135 == half_angle(  0, 270)

        assert  -45 == half_angle( 90,   0)
        assert  -90 == half_angle(180,   0)
        assert -135 == half_angle(270,   0)

        assert -135 == half_angle( -90%360, 0) # 270 -> 0
        assert  -90 == half_angle(-180%360, 0) # 180 -> 0
        assert  -45 == half_angle(-270%360, 0) #  90 -> 0

        assert -45 == half_angle(180,  90)
        assert   0 == half_angle(180, 180)
        assert  45 == half_angle(180, 270)

        # NOTE: It is chosen to be consistently this way ...
        assert -90 == half_angle(180,   0)
        assert  90 == half_angle(180,   0, fwd = True)

        assert  90 == half_angle(  0, 180)
        assert  95 == half_angle(  0, 190)
        assert  90 == half_angle( 10, 190)

    def test_get_coords(self):
        assert np.allclose((1, 0), get_coords((0,0),    0, 1))
        assert np.allclose((0, 1), get_coords((0,0),   90, 1)) 
        assert np.allclose((-1, 0), get_coords((0,0), 180, 1))
        assert np.allclose((0, -1), get_coords((0,0), 270, 1)) 

        assert np.allclose((2, 1), get_coords((1,1),   0, 1))
        assert np.allclose((1, 2), get_coords((1,1),  90, 1)) 
        assert np.allclose((0, 0), get_coords((1,0), 180, 1))
        assert np.allclose((1, -1), get_coords((1,0), 270, 1)) 

    def test_get_a1(self):
        # Missing previous angle
        assert 150 == get_a1(None, 0)
        assert 180 == get_a1(None, 30)
        assert 240 == get_a1(None, 90)
        assert 330 == get_a1(None, 180)
        assert  60 == get_a1(None, 270)

        # Given angle = 0
        assert -90+180 == get_a1(180, 0) # 180 points up!
        assert -89+180 == get_a1(182, 0)
        assert -45+180 == get_a1(270, 0)
        assert  -1+180 == get_a1(358, 0)
        assert     180 == get_a1(  0, 0)
        assert   1+180 == get_a1(  2, 0)
        assert  45+180 == get_a1( 90, 0)
        assert  89+180 == get_a1(178, 0)

        # Given angle = 90
        assert   1+180 == get_a1(272, 90)
        assert  45+180 == get_a1(  0, 90)
        assert  90+180 == get_a1( 90, 90)
        assert 135+180 == get_a1(180, 90)
        assert 179+180 == get_a1(268, 90)
        assert     180 == get_a1(270, 90)

    def test_get_a2(self):
        # Missing previous angle
        assert  30 == get_a2(  0, None)
        assert  20 == get_a2(350, None)

        # Given angle = 0
        assert   0 == get_a2(0,   0)
        assert   1 == get_a2(0,   2)
        assert  45 == get_a2(0,  90)
        assert  89 == get_a2(0, 178)
        assert  90 == get_a2(0, 180) # 180 points up!
        assert 271 == get_a2(0, 182)
        assert 315 == get_a2(0, 270)
        assert 359 == get_a2(0, 358)

        # Given angle = 90
        assert  45 == get_a2(90,   0)
        assert  90 == get_a2(90,  90)
        assert 135 == get_a2(90, 180)
        assert 179 == get_a2(90, 268)
        assert 180 == get_a2(90, 270)
        assert   1 == get_a2(90, 272)

    def test_get_a3(self):
        assert 150 == get_a3(None, 0)
        assert  91+180 == get_a3(  2, 180)
        assert  90+180 == get_a3(  0, 180)
        assert  89 == get_a3(358, 180)
        assert  45 == get_a3(270, 180)
        assert   0 == get_a3(180, 180)
        assert -45+360 == get_a3( 90, 180)

    def test_get_a4(self):
        assert 210 == get_a4(180, None)
        assert -89+180 == get_a4(180,   2)
        assert  90+180 == get_a4(180,   0)
        assert  89+180 == get_a4(180, 358)
        assert  45+180 == get_a4(180, 270)
        assert     180 == get_a4(180, 180)
        assert -45+180 == get_a4(180,  90)

