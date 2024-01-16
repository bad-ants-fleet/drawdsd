#
# Produce an SVG file from data.
#
import logging
log = logging.getLogger(__name__)

import numpy as np
import drawsvg as draw
from colorsys import hls_to_rgb
from .components import scale, agl
from scipy.optimize import fsolve


dfont = 'bold' # TODO: need some interfacd for plotting parameters.
ffamily = 'Courier'

# Layout
def dom_path(dcolor, sw = None, da = None):
    if sw is None:
        sw = 10 
    return draw.Path(stroke = dcolor, 
                     fill_opacity = 0, 
                     stroke_width = sw, 
                     stroke_dasharray = da,
                     stroke_linejoin = "round")

def dom_txt(dname, x, y, ux, uy):
    x = x - uy*scale
    y = y + ux*scale
    return draw.Text(dname, 14, x = round(x), y = round(y), 
                     font_family = ffamily,
                     font_weight = dfont, 
                     dominant_baseline='middle',
                     text_anchor='middle')

def h_bonds(x, y, ux, uy):
    p = draw.Path(stroke = 'black',
                  fill_opacity = 0, 
                  stroke_width = 2)
    x0 = x + uy * scale/2
    y0 = y - ux * scale/2
    xn = x + uy * scale*1.2
    yn = y - ux * scale*1.2
    p.M(round(x0), round(y0)).L(round(xn), round(yn))
    return p

def nuc_circle(x, y, c, ux, uy, off = scale/8):
    x = x + uy*off
    y = y - ux*off
    return draw.Circle(x, y, scale/3, fill = c, stroke_width = 0, stroke = 'black')

def nuc_txt(char, x, y, ux, uy, off = scale/8):
    x = x + uy*off
    y = y - ux*off
    return draw.Text(char, round(scale*2/3), round(x), round(y), 
                     font_family = ffamily,
                     font_weight = dfont, 
                     dominant_baseline='middle',
                     text_anchor='middle')

#def draw_5prime(p, xs, ys, ux, uy):
#    a = 1.5
#    b = 0.5
#    c = (a+b)/2
#    x50, y50 = xs  + ux*scale*a, ys + uy*scale*a
#    x51, y51 = xs  + ux*scale*c, ys + uy*scale*c
#    x52, y52 = x51 - uy*scale, y51  + ux*scale
#    x53, y53 = xs  + ux*scale*b, ys + uy*scale*b
#    return p.M(round(x50), round(y50)).Q(x52, y52, x53, y53)  # Draw a curve to (70, -20)

def draw_5prime(p, xs, ys, ux, uy):
    a, b, c = 1, 0.5, 0.5
    x50, y50 = xs  + ux*scale*a, ys  + uy*scale*a
    x51, y51 = x50 - uy*scale*b, y50 + ux*scale*b
    x52, y52 = x51 - ux*scale*c, y51 - uy*scale*c
    x53, y53 = x52 + uy*scale*b, y52 - ux*scale*b
    return p.M(round(x50), round(y50)).L(
               round(x51), round(y51)).L(
               round(x52), round(y52)).L(
               round(x53), round(y53))

def draw_3prime(p, xn, yn, ux, uy):
    a, b, c = 1.5, 1, 1
    rx, ry = vec_rot(ux, uy, 180)
    x30, y30 = xn + rx*scale*a, yn + ry*scale*a
    x31, y31 = x30 + ux*scale*b, y30 + uy*scale*b
    rx, ry = vec_rot(ux, uy, 90)
    x32, y32 = x30 + rx*scale*c, y30 + ry*scale*c
    return p.L(round(x30), round(y30)).L(
               round(x31), round(y31)).L(
               round(x32), round(y32))

def get_radius_2(s, c, points):
    def frad(x): # minor arc
        (r, a) = x
        return [s/points - (2 * r * np.sin((a/points)/(2*r))),
                np.sin(a/(2*r)) - c/(2*r)]
    # Initial guess: assume a+c = U
    [r, a] = fsolve(frad, [(s+c)/(2*np.pi), s])
    assert 2*r > c
    if r*np.pi > a: # minor
        phi = agl(np.rad2deg(a/r))
        psi = phi
        off = -(180-phi)/2
    else: # mayor
        phi = 360 - agl(np.rad2deg(a/r))
        psi = 360 - phi
        off = (180-phi)/2
    return r, phi, psi, off


def get_radius(a, c):
    assert a > c
    def frad(r): # minor arc
        return np.sin(a/(2*r)) - c/(2*r)
    # Initial guess: assume a+c = U
    [r] = fsolve(frad, (a+c)/(2*np.pi))
    assert 2*r > c

    if r*np.pi > a: # minor
        phi = agl(np.rad2deg(a/r))
        psi = phi
        off = -(180-phi)/2
    else: # mayor
        phi = 360 - agl(np.rad2deg(a/r))
        psi = 360 - phi
        off = (180-phi)/2
    return r, phi, psi, off

#def get_offset(detour, dist_0N):
#    # The angle assuming a full circle
#    ang = 180 - 360/((detour + dist_0N)/scale)
#    print(f'{ang=}')
#    mang = ((180 - ang)*dist_0N/scale)
#    print(f'{mang=}')
#    return 180 - mang

def get_detour_circ(p0, pN, dtarget):
    dist_0N = distance(p0, pN) # chordlength
    pts = [p0]
    apN = True
    if a_gt_b((dtarget - 2*scale), dist_0N):
        # Chop two points and put them separately
        tux, tuy = unitvec(p0, pN, dist_0N)
        rux, ruy = vec_rot(tux, tuy, 90)
        pl0 = end_point(p0, (rux, ruy), scale)
        plN = end_point(pN, (rux, ruy), scale)
        ddetour = dtarget - 2*scale
        pts.append(pl0)
    else:
        pl0 = p0
        plN = pN
        ddetour = dtarget
        apN = False

    print('--')
    c = dist_0N # cordlength
    a = ddetour  # arclength
    points = round(a/scale-0.5)
    print(f'{c = }, {a = }')

    assert not a_eq_b(a, c)

    r, phi, psi, off = get_radius(a, c)
    print(f'{r=} {phi=} {psi=} {off=}')
    assert a_eq_b(np.deg2rad(psi)*r, a)
    r, phi, psi, off = get_radius_2(a, c, points)
    print(np.deg2rad(psi)*r, a)
    print(f'{r=} {phi=} {psi=} {off=}')

    # The center calculated from p0
    tux, tuy = unitvec(pl0, plN)
    rux, ruy = vec_rot(tux, tuy, off)
    cc = (end_point(pl0, (rux, ruy),  r))

    # The center calculated from pN
    tux, tuy = unitvec(plN, pl0)
    rux, ruy = vec_rot(tux, tuy, -off)
    cc2 = (end_point(plN, (rux, ruy),  r))
    assert np.allclose(cc, cc2)

    sa = a / points
    sc = 2*r*np.sin(sa/(2*r))
    print(sa, sc, points*sc, a)

    ang = psi/points
    for _ in range(points):
        #if _ == 0:
        #    tux, tuy = unitvec(pl1, cc1)
        #    rux, ruy = vec_rot(tux, tuy, ang)
        #    pts.append(end_point(pts[-1], (rux, ruy),  a))
        #else:
        tux, tuy = unitvec(cc, pts[-1])
        rux, ruy = vec_rot(tux, tuy, -ang)
        pts.append(end_point(cc, (rux, ruy),  r))
    if apN:
        pts.append(pN)
    return pts#(*cc1, r), (*cc2, r)


    ## Target segment length
    ##print(f'{ddetour = } {points = }')
    #a = ddetour/points
    ##print(f'a {a=}')
    ##tphi = a/r
    ##print(f' {tphi=}')
    ##a = abs(2*r*np.sin(tphi/2))
    ##print(f'b {a=}')
    #print(distance(p0, pN))
    #print(distance(pl1, plN))
    #print(distance(pl1, pts[-1]), distance(pl1, plN))
    #print(distance(plN, pts[-1]), distance(plN, plN))
    #if True:
    #    #print(pl1)
    #    #print(pts[-1])
    #    #print(plN)
    #    # Center piece
    #    tux, tuy = unitvec(pts[-1], pts[-2])
    #    rux, ruy = vec_rot(tux, tuy, ang/2 + 90)
    #    pts.append(end_point(pts[-1], (rux, ruy),  dist_0N))
    #    # Second half 
    #    for _ in range(points):
    #        if _ == 0:
    #            tux, tuy = unitvec(pts[-1], pts[-2])
    #            rux, ruy = vec_rot(tux, tuy, 90+ang/2)
    #            pts.append(end_point(pts[-1], (rux, ruy),  a))
    #        else:
    #            tux, tuy = unitvec(pts[-1], pts[-2])
    #            rux, ruy = vec_rot(tux, tuy, ang)
    #            pts.append(end_point(pts[-1], (rux, ruy),  a))
    #else:

    return (0,0)


def get_detour(p0, pN, dtarget, circ = False):
    dist_0N = distance(p0, pN)
    
    pts = [p0]
    apN = True
    if a_ge_b((dtarget - 2*scale), dist_0N):
        # Chop two points and put them separately
        tux, tuy = unitvec(p0, pN, dist_0N)
        rux, ruy = vec_rot(tux, tuy, 90)
        pl1 = end_point(p0, (rux, ruy), scale)
        plN = end_point(pN, (rux, ruy), scale)
        ddetour = dtarget - 2*scale - dist_0N
        pts.append(pl1)
    else:
        pl1 = p0
        plN = pN
        ddetour = dtarget - dist_0N
        apN = False

    points = round(ddetour/scale/2-0.5) if circ else 1
    a = ddetour/2/points
    ang = 180 - (360/(points*2))
    for _ in range(points):
        if _ == 0:
            tux, tuy = unitvec(pts[-1], plN)
            rux, ruy = vec_rot(tux, tuy, 90+ang/2)
            pts.append(end_point(pts[-1], (rux, ruy),  a))
        else:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, ang)
            pts.append(end_point(pts[-1], (rux, ruy),  a))
    # Center piece
    tux, tuy = unitvec(pts[-1], pts[-2])
    rux, ruy = vec_rot(tux, tuy, ang/2 + 90)
    pts.append(end_point(pts[-1], (rux, ruy),  dist_0N))
    # Second half 
    for _ in range(points):
        if _ == 0:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, 90+ang/2)
            pts.append(end_point(pts[-1], (rux, ruy),  a))
        else:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, ang)
            pts.append(end_point(pts[-1], (rux, ruy),  a))
    if apN:
        pts.append(pN)
    return pts

# Helpers
def a_lt_b(a, b):
    return a < b and not np.isclose(a, b)

def a_eq_b(a, b, eps=1e-5):
    return np.isclose(a, b, rtol=eps)

def a_gt_b(a, b):
    return a > b and not np.isclose(a, b)

def a_ge_b(a, b):
    return a > b or np.isclose(a, b)


def distance(p1, p2):
    """Returns distance between two points."""
    return np.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

def unitvec(p1, p2, dist = None):
    """Returns unit vector from p1 to p2."""
    if dist is None:
        dist = distance(p1, p2) 
    return (p2[0]-p1[0])/dist, (p2[1]-p1[1])/dist

def vec_rot(x, y, alpha):
    """Rotates vector by angle alpha. """
    alpha = np.radians(alpha)
    rx = x * np.cos(alpha) - y * np.sin(alpha)
    ry = x * np.sin(alpha) + y * np.cos(alpha)
    return rx, ry

def end_point(p0, v, l):
    """Returns the end point given p0 vector. """
    return p0[0] + v[0]*l, p0[1] + v[1]*l

# Interface 
def get_drawing(svgC):
    # Initialize the SVG image & background
    dimx, dimy, minx, miny = estimate_dimensions(svgC)

    # Add a 50px frame around the molecule
    dimx += 100
    dimy += 100
    minx -= 50
    miny -= 50
    
    svg = draw.Drawing(dimx, dimy, 
                       origin = (minx, miny), displayInline = False)
    svg.append(draw.Rectangle(minx, miny, 
                              dimx, dimy, fill='white'))
 
    svg.extend(svgC)
    return svg

def estimate_dimensions(svgC):
    minx, maxx, miny, maxy = 0, 0, 0, 0
    for obj in svgC:
        if 'x' in obj.args:
            if obj.args['x'] > maxx:
                maxx =  obj.args['x']
            if obj.args['x'] < minx:
                minx =  obj.args['x']
            if obj.args['y'] > maxy:
                maxy =  obj.args['y']
            if obj.args['y'] < miny:
                miny =  obj.args['y']
        else:
            try:
                for x, y in (x[1:].split(',') for x in obj.args['d'].split()):
                    if int(x) > maxx:
                        maxx = int(x)
                    if int(x) < minx:
                        minx = int(x)
                    if int(y) > maxy:
                        maxy = int(y)
                    if int(y) < miny:
                        miny = int(y)
            except KeyError:
                pass
            except ValueError:
                pass
    dimx = maxx - minx
    dimy = maxy - miny
    return dimx, dimy, minx, miny

def get_rgb_palette(num):
    num = int(360/(num+1))
    palette = [hls_to_rgb(angle/360, .6, .8) for angle in range(0, 360, num)]
    return [(int(r*255), int(g*255), int(b*255)) for (r,g,b) in palette]

def draw_stem(i1, i2, i3, i4, p15, p23, p35, p43, color, angle, length, nameT, nameB):
    circles = []
    bground = draw.Path(stroke_width=0, fill='black', fill_opacity=.3)
    dompath = dom_path(color)

    (x1, y1) = i1
    (x2, y2) = i2
    (x3, y3) = i3
    (x4, y4) = i4

    bground.M(round(x1), round(y1)).L(
              round(x2), round(y2)).L(
              round(x3), round(y3)).L(
              round(x4), round(y4))

    tx, ty = (x1+x2)/2, (y1+y2)/2 
    ux, uy = unitvec((x1, y1), (x2, y2))
    dt1 = dom_txt(nameT, tx, ty, ux, uy)
    cx, cy = x1 + ux*scale/2, y1 + uy*scale/2
    for _ in range(length):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt('N', cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
 
    if p15:
        dompath = draw_5prime(dompath, x1, y1, ux, uy)
    else:
        dompath.M(round(x1), round(y1))

    if p23:
        dompath = draw_3prime(dompath, x2, y2, ux, uy)
    else:
        dompath.L(round(x2), round(y2))

    tx, ty = (x3+x4)/2, (y3+y4)/2 
    ux, uy = unitvec((x3, y3), (x4, y4))
    dt2 = dom_txt(nameB, tx, ty, ux, uy)
    cx, cy = x3 + ux*scale/2, y3 + uy*scale/2
    for _ in range(length):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt('N', cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
 
    if p35:
        dompath = draw_5prime(dompath, x3, y3, ux, uy)
    else:
        dompath.M(round(x3), round(y3))

    if p43:
        dompath = draw_3prime(dompath, x4, y4, ux, uy)
    else:
        dompath.L(round(x4), round(y4))

    return bground, dompath, dt1, dt2, *circles

def draw_tentacles(dns, dls, dcs, dzs, dos, dms):
    """
    """
    # names, lengths, colors, startpoint, endpoint
    for ns, ls, cs, (x0, y0), (x1, y1), es in zip(dns, dls, dcs, dzs, dos, dms):
        slength = scale * sum(ls) # total length of concatenated backbones
        linelen = distance((x0, y0), (x1, y1))
        if a_eq_b(linelen, 0):
            continue # Necessary 
        #yield draw.Circle(0, 0, 285.14, stroke_width = 1, fill = None),
        #yield draw.Circle(0, 0, 2, stroke_width = 1, fill = 'white'),
        #tp = draw.Path(stroke='white')
        #tp0 = (0, 0)
        #tp1 = (0, 285.15)
        #tp2 = vec_rot(*tp1, 80.374)
        #tp.M(*tp1).L(*tp0).L(*tp2)
        #yield tp,

        # Append a fake domain if the distance between points is 
        # longer than the length of all domains.
        if a_lt_b(slength, linelen):
            ns.append('')
            ls.append((linelen-slength)/scale)
            cs.append('black')
            es.append(False)
            slength = scale * sum(ls)
            assert a_eq_b(slength, linelen)

        if a_eq_b(slength, linelen):
            ux, uy = unitvec((x0, y0), (x1, y1), linelen)
            xs, ys = x0, y0 # start coordinates
            circles = []
            for (dname, dlen, dcolor, end) in zip(ns, ls, cs, es):
                # Initialize path
                sw, da = (None, None) if dname else (2, 5) # fake domain
                p = dom_path(dcolor, sw, da)
                sdlen = scale * dlen
                xn, yn = xs + ux*sdlen, ys + uy*sdlen
                if end == 'p5':
                    p = draw_5prime(p, xs, ys, ux, uy)
                    #x5, y5 = xs  + ux*scale/1.5, ys + uy*scale/1.5
                    #yield draw.Circle(x5, y5, scale/1.5, stroke_width = 0, fill = dcolor),
                    p.L(round(xn), round(yn))
                elif end == 'p3':
                    p.M(round(xs), round(ys))
                    p = draw_3prime(p, xn, yn, ux, uy)
                else:
                    p.M(round(xs), round(ys))
                    p.L(round(xn), round(yn))
                # Draw nametag
                tx, ty = (xs+xn)/2, (ys+yn)/2 # point in the middle on straight line
                dt = dom_txt(dname, tx, ty, ux, uy)
                if dname:
                    cx, cy = xs + ux*scale/2, ys + uy*scale/2
                    for _ in range(dlen):
                        circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                        circles.append(nuc_txt('N', cx, cy, ux, uy))
                        cx, cy = cx + ux*scale, cy + uy*scale
                xs, ys = xn, yn
                yield p, dt
                yield circles

        else: # not enough space, let's draw a rectangle!
            assert slength > linelen
            detour = get_detour_circ((x0, y0), (x1, y1), slength)
            #cp1, cp2 = get_detour_circ((x0, y0), (x1, y1), slength)
            #yield [draw.Circle(*cp1, stroke_width = 1, fill = 'black')]
            #yield [draw.Circle(*cp2, stroke_width = 1, fill = 'white')]
            #yield [draw.Circle(x0, y0, 5, stroke_width = 1, fill = 'black')]
            #yield [draw.Circle(x1, y1, 5, stroke_width = 1, fill = 'black')]

            #    yield [draw.Circle(*p, scale/2, stroke_width = 1, stroke = 'black', fill = 'black')]
            circles = []
            (xs, ys), di = detour[0], 1
            off = scale/2
            for (dname, dlen, dcolor, end) in zip(ns, ls, cs, es):
                rsdlen = scale * dlen # remaining scaled domain length
                tlabel, tl = scale * dlen / 2, True
                p = dom_path(dcolor)
                p.M(round(xs), round(ys))
                while True:
                    xn, yn = detour[di]
                    seglen = distance((xs, ys), (xn, yn))
                    ux, uy = unitvec((xs, ys), (xn, yn), seglen)
                    # Domain label
                    if tl and (np.isclose(seglen, tlabel) or seglen > tlabel):
                        tx, ty = xs + ux*tlabel, ys + uy*tlabel
                        dt = dom_txt(dname, tx, ty, ux, uy)
                        tl = False
                    # Circle initialization
                    cx, cy = xs + ux*off, ys + uy*off
                    if a_gt_b(distance((xs, ys), (cx, cy)), seglen): # a >= b
                        off = distance((xs, ys), (cx, cy)) - seglen
                    # Finish the path if segment length >= remaining domain length
                    if a_eq_b(rsdlen, seglen) or seglen > rsdlen:
                        xn, yn = xs + ux*rsdlen, ys + uy*rsdlen
                        # Draw domain
                        if end == 'p5':
                            p = draw_5prime(p, xs, ys, ux, uy)
                            p.L(round(xn), round(yn))
                            end = False
                        elif end == 'p3':
                            p.M(round(xs), round(ys))
                            p = draw_3prime(p, xn, yn, ux, uy)
                            end = False
                        else:
                            p.L(round(xs), round(ys))
                            p.L(round(xn), round(yn))
                        # Draw nucleotides
                        for _ in range(dlen):
                            circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                            circles.append(nuc_txt('N', cx, cy, ux, uy),)
                            cx, cy = cx + ux*scale, cy + uy*scale
                            if a_gt_b(distance((xs, ys), (cx, cy)), rsdlen): # a >= b
                                off = scale/2
                                break
                        if a_eq_b(rsdlen, seglen):
                            di+=1
                        xs, ys = xn, yn
                        break
                    # Draw domain
                    if end == 'p5':
                        p = draw_5prime(p, xs, ys, ux, uy)
                        p.L(round(xn), round(yn))
                        end = False
                    else:
                        p.L(round(xs), round(ys))
                        p.L(round(xn), round(yn))
                    # Draw nucleotides
                    for _ in range(dlen):
                        circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                        circles.append(nuc_txt('N', cx, cy, ux, uy))
                        cx, cy = cx + ux*scale, cy + uy*scale
                        if a_gt_b(distance((xs, ys), (cx, cy)), seglen): # a >= b
                            off = distance((xs, ys), (cx, cy)) - seglen
                            break
                    rsdlen -= seglen
                    if tl:
                        tlabel -= seglen
                    di += 1 
                    xs, ys = xn, yn
                yield p, dt
            yield circles

