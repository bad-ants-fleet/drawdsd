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

# Global parameters
prec = 4
dfont = 'bold' # TODO: need some interfacd for plotting parameters.
ffamily = 'Helvetica'

# Layout
def dom_path(dcolor, sw = None, da = None):
    if sw is None:
        sw = 9
    p = draw.Path(stroke = dcolor, 
                     fill_opacity = 0, 
                     stroke_width = sw, 
                     stroke_dasharray = da,
                     stroke_linecap='butt',
                     stroke_linejoin = "round")
    p.ddsdinfo = 'backbone'
    return p


def dom_txt(dname, x, y, ux, uy):
    x = x + uy*scale
    y = y - ux*scale
    return draw.Text(dname, 14, x = round(x, prec), y = round(y, prec), 
                     font_family = ffamily,
                     font_weight = dfont, 
                     dominant_baseline='middle',
                     text_anchor='middle')

def h_bonds(x, y, ux, uy):
    p = draw.Path(stroke = 'black',
                  fill_opacity = 0, 
                  stroke_width = 2)
    x0 = x - uy * scale/2
    y0 = y + ux * scale/2
    xn = x - uy * scale*1.2
    yn = y + ux * scale*1.2
    p.M(round(x0, prec), round(y0, prec)).L(round(xn, prec), round(yn, prec))
    p.ddsdinfo = 'base-pair'
    return p

def nuc_circle(x, y, c, ux, uy, off = scale/8):
    x = x - uy*off
    y = y + ux*off
    return draw.Circle(round(x, prec), round(y, prec), round(scale/3, prec), fill = c, stroke_width = 0, stroke = 'black')

def nuc_txt(char, x, y, ux, uy, off = scale/8):
    x = x - uy*off
    y = y + ux*off
    return draw.Text(char, round(scale*2/3, prec), round(x, prec), round(y, prec), 
                     font_family = ffamily,
                     font_weight = dfont, 
                     dominant_baseline='middle',
                     text_anchor='middle')

#def draw_5prime(p, xs, ys, ux, uy):
#    a, b = 1.5, 0.5 
#    c = (a+b)/2
#    x50, y50 = xs  + ux*scale*a, ys + uy*scale*a
#    x51, y51 = xs  + ux*scale*c, ys + uy*scale*c
#    x52, y52 = x51 - uy*scale, y51  + ux*scale
#    x53, y53 = xs  + ux*scale*b, ys + uy*scale*b
#    return p.M(round(x50), round(y50)).Q(x52, y52, x53, y53)  # Draw a curve to (70, -20)

def draw_5prime(p, xs, ys, ux, uy):
    a, b, c = 1, .2, 1
    #x50, y50 = xs + uy*scale*a, ys - ux*scale*a
    x51, y51 = xs + ux*scale*b, ys + uy*scale*b
    #x52, y52 = x51 - ux*scale*c, y51 - uy*scale*c
    #x53, y53 = x52 + uy*scale*b, y52 - ux*scale*b
    #p.M(round(x50), round(y50))
    return p.M(round(x51, prec), round(y51, prec))
#.L(
            #   round(x52), round(y52))
            #   #.L( round(x53), round(y53))

def draw_3prime(p, xn, yn, ux, uy):
    a, b, c, d = 1.1, 0.7, 1, 1.1,
    rx, ry = vec_rot(ux, uy, 180)
    x30, y30 = xn + rx*scale*a, yn + ry*scale*a
    x31, y31 = x30 + ux*scale*b, y30 + uy*scale*b
    x33, y33 = xn + rx*scale*d, yn + ry*scale*d
    rx, ry = vec_rot(ux, uy, -90)
    x32, y32 = x33 + rx*scale*c, y33 + ry*scale*c
    return p.L(round(x30, prec), round(y30, prec)).L(
               round(x31, prec), round(y31, prec)).L(
               round(x32, prec), round(y32, prec))

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
        rux, ruy = vec_rot(tux, tuy, -90)
        pl0 = end_point(p0, (rux, ruy), scale)
        plN = end_point(pN, (rux, ruy), scale)
        ddetour = dtarget - 2*scale
        pts.append(pl0)
    else:
        pl0 = p0
        plN = pN
        ddetour = dtarget

    c = dist_0N # cordlength
    a = ddetour # arclength
    points = 10 *round(a/scale)
    assert not a_eq_b(a, c)

    r, phi, psi, off = get_radius(a, c)
    assert a_eq_b(np.deg2rad(psi)*r, a)

    # Estimate error from lengths
    sa = a / points
    sc = 2*r*np.sin(sa/(2*r))
    err = abs(points*sc - a)

    # The center calculated from p0
    tux, tuy = unitvec(pl0, plN)
    rux, ruy = vec_rot(tux, tuy, off)
    cc = (end_point(pl0, (rux, ruy),  r))

    # The center calculated from pN
    tux, tuy = unitvec(plN, pl0)
    rux, ruy = vec_rot(tux, tuy, -off)
    cc2 = (end_point(plN, (rux, ruy),  r))
    assert np.allclose(cc, cc2)

    ran = 180 - (360 - psi)
    ang = 180 - psi/points
    for i in range(points):
        d = sc + err/(points)
        if i == 0:
            tux, tuy = unitvec(pl0, plN)
            rux, ruy = vec_rot(tux, tuy, -(ran/2 + ang/2))
            pts.append(end_point(pts[-1], (rux, ruy),  d))
        else:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy,- ang)
            pts.append(end_point(pts[-1], (rux, ruy),  d))

    if pts[-1] != plN:
        ed = distance(pts[-1], plN)
        tux, tuy = unitvec(pts[-1], plN)
        newp = [p0]
        for p in pts:
            p = end_point(p, (tux, tuy),  ed/2)
            newp.append(p)
        pts = newp
    
    pts.append(pN) # always append target coordinate

    return pts

def get_detour(p0, pN, dtarget, circ = False):
    dist_0N = distance(p0, pN)
    
    pts = [p0]
    apN = True
    if a_ge_b((dtarget - 2*scale), dist_0N):
        # Chop two points and put them separately
        tux, tuy = unitvec(p0, pN, dist_0N)
        rux, ruy = vec_rot(tux, tuy, -90)
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
            rux, ruy = vec_rot(tux, tuy, -(90+ang/2))
            pts.append(end_point(pts[-1], (rux, ruy),  a))
        else:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, -ang)
            pts.append(end_point(pts[-1], (rux, ruy),  a))
    # Center piece
    tux, tuy = unitvec(pts[-1], pts[-2])
    rux, ruy = vec_rot(tux, tuy, -(ang/2 + 90))
    pts.append(end_point(pts[-1], (rux, ruy),  dist_0N))
    # Second half 
    for _ in range(points):
        if _ == 0:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, -(90+ang/2))
            pts.append(end_point(pts[-1], (rux, ruy),  a))
        else:
            tux, tuy = unitvec(pts[-1], pts[-2])
            rux, ruy = vec_rot(tux, tuy, -ang)
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
    svg.append(draw.Rectangle(minx+20, miny+20, 
                              dimx-40, dimy-40, stroke_width=5, stroke = 'black', fill='none'))
    
    bg_layer = []
    backbone = []
    nuc_layer = []
    txt_layer = []
    for x in svgC:
        if x.TAG_NAME == 'text':
            txt_layer.append(x)
        elif x.TAG_NAME == 'circle':
            nuc_layer.append(x)
        elif x.TAG_NAME == 'path':
            if x.ddsdinfo == 'bground':
                bg_layer.append(x)
            elif x.ddsdinfo == 'base-pair':
                nuc_layer.append(x)
            elif x.ddsdinfo == 'backbone':
                backbone.append(x)
            else:
                print(f'Unknown path object: {x = }')
        else:
            print(f'Unknown object TAG: {x = }')

    svg.extend(bg_layer)

    bbn = []
    for b in backbone:
        p = draw.Path(**b.args)
        p.args['d'] = truncate_domains(p.args['d'], True, True)
        bbn.append(p)

    def get_ends(seg):
        ms = seg.args['d'].split()
        assert ms[0][0] == 'M'
        p0 = tuple(ms[0][1:].split(','))
        assert ms[-1][0] == 'L'
        pN = tuple(ms[-1][1:].split(','))
        return p0, pN

    sp0s = {}
    spNs = {}
    segd = {}
    for s in backbone:
        p0, pN = get_ends(s)
        segd[(p0, pN)] = s
        sp0s[p0] = sp0s.get(p0, []) + [s]
        spNs[pN] = sp0s.get(pN, []) + [s]

    paths = []
    while len(backbone):
        seg = backbone.pop(0)
        path = [seg]
        p0, pN = get_ends(seg)
        seen = []
        while len(backbone):
            nseg = backbone.pop(0)
            np0, npN = get_ends(nseg)
            if npN == p0:
                path.insert(0, nseg)
                p0 = np0
                seen = []
            elif np0 == pN:
                path.append(nseg)
                pN = npN
                seen = []
            else:
                backbone.append(nseg)
                seen.append(nseg)
            if len(seen) == len(backbone):
                break
        paths.append(path)

    for path in paths:
        p = draw.Path(**(path[0]).args)
        p.args['stroke'] = 'gray'
        p.args['stroke-width'] = 16
        p.args['stroke-linecap'] = 'butt'
        d = p.args['d']
        for seg in path[1:]:
            x = seg.args['d'].split()
            y = ' '.join(x[1:])
            p.args['d'] += ' '+y
        #p.args['d'] = truncate_domains(p.args['d'], True, True)
        svg.append(p)

    svg.extend(bbn)
    svg.extend(nuc_layer)
    svg.extend(txt_layer)
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
    palette = [hls_to_rgb(angle/360, .7, .8) for angle in range(0, 360, num)]
    return [(int(r*255), int(g*255), int(b*255)) for (r,g,b) in palette]

def draw_stem(i1, i2, i3, i4, p15, p23, p35, p43, color, angle, length, seqT, seqB, nameT, nameB):
    circles = []
    dompathT = dom_path(color)
    tseqT = list(seqT)
    tseqB = list(seqB)

    (x1, y1) = i1
    (x2, y2) = i2
    (x3, y3) = i3
    (x4, y4) = i4

    bground = draw.Lines(x1, y1, x2, y2, x3, y3, x4, y4,
                stroke='none', fill='lightgray', close='true')
    bground.ddsdinfo = 'bground'

    tx, ty = (x1+x2)/2, (y1+y2)/2 
    ux, uy = unitvec((x1, y1), (x2, y2))
    dt1 = dom_txt(nameT, tx, ty, ux, uy)
    cx, cy = x1 + ux*scale/2, y1 + uy*scale/2
    while len(tseqT):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt(tseqT.pop(0), cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
 
    #if p15:
    #    ux, uy = unitvec((x1, y1), (x2, y2))
    #    dompathT = draw_5prime(dompathT, x1, y1, ux, uy)
    #else:
    dompathT.M(round(x1, prec), round(y1, prec))

    if p23:
        dompathT = draw_3prime(dompathT, x2, y2, ux, uy)
    else:
        dompathT.L(round(x2, prec), round(y2, prec))

    dompathB = dom_path(color)
    tx, ty = (x3+x4)/2, (y3+y4)/2 
    ux, uy = unitvec((x3, y3), (x4, y4))
    dt2 = dom_txt(nameB, tx, ty, ux, uy)
    cx, cy = x3 + ux*scale/2, y3 + uy*scale/2
    while len(tseqB):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt(tseqB.pop(0), cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
 
    #if p35:
    #    ux, uy = unitvec((x1, y1), (x2, y2))
    #    dompathB = draw_5prime(dompathB, x3, y3, ux, uy)
    #else:
    dompathB.M(round(x3, prec), round(y3, prec))

    if p43:
        dompathB = draw_3prime(dompathB, x4, y4, ux, uy)
    else:
        dompathB.L(round(x4, prec), round(y4, prec))

    return bground, dompathT, dompathB, dt1, dt2, *circles


def truncate_domains(d, p5, p3):
    x = [[c[0], *c[1:].split(',')] for c in d.split()]
    if p5:
        x1 = float(x[0][1])
        y1 = float(x[0][2])
        x2 = float(x[1][1])
        y2 = float(x[1][2])
        ux, uy = unitvec((x1, y1), (x2, y2))
        myd = min(distance((x1,y1),(x2,y2)), scale/6)
        ps = end_point((x1, y1), (ux, uy), myd)
        x[0][1], x[0][2] = str(round(ps[0], prec)), str(round(ps[1], prec))
    if p3:
        xs = float(x[-2][1])
        ys = float(x[-2][2])
        xn = float(x[-1][1])
        yn = float(x[-1][2])
        ux, uy = unitvec((xn, yn), (xs, ys))
        myd = min(distance((x1,y1),(x2,y2)), scale/6)
        ps = end_point((xn, yn), (ux, uy), myd)
        x[-1][1], x[-1][2] = str(round(ps[0], prec)), str(round(ps[1], prec))
    y = ' '.join([c[0]+','.join(c[1:]) for c in x])
    return y

def draw_tentacles(dns, dls, dcs, dzs, dos, dss, dms):
    """
    """
    # names, lengths, colors, startpoint, endpoint
    for ns, ls, cs, (x0, y0), (x1, y1), ss, es in zip(dns, dls, dcs, dzs, dos, dss, dms):
        slength = scale * sum(ls) # total length of concatenated backbones
        linelen = distance((x0, y0), (x1, y1))
        if a_eq_b(linelen, 0):
            assert slength == 0
            continue # No tentacle.
        # Append a fake domain if the distance between points is 
        # longer than the length of all domains.
        if a_lt_b(slength, linelen):
            ns.append('')
            ls.append((linelen-slength)/scale)
            cs.append('black')
            es.append(False)
            ss.append('')
            slength = scale * sum(ls)
            assert a_eq_b(slength, linelen)
        if a_eq_b(slength, linelen):
            ux, uy = unitvec((x0, y0), (x1, y1), linelen)
            xs, ys = x0, y0 # start coordinates
            circles = []
            for (dname, dlen, dcolor, seq, end) in zip(ns, ls, cs, ss, es):
                # Initialize path
                tseq = list(seq)
                sw, da = (None, None) if dname else (2, 5) # fake domain
                p = dom_path(dcolor, sw, da)
                sdlen = scale * dlen
                xn, yn = xs + ux*sdlen, ys + uy*sdlen
                #if end == 'p5':
                #    p = draw_5prime(p, xs, ys, ux, uy)
                #    #x5, y5 = xs  + ux*scale/1.5, ys + uy*scale/1.5
                #    #yield draw.Circle(x5, y5, scale/1.5, stroke_width = 0, fill = dcolor),
                #    p.L(round(xn, prec), round(yn, prec))
                #
                if end == 'p3':
                    p.M(round(xs, prec), round(ys, prec))
                    p = draw_3prime(p, xn, yn, ux, uy)
                else:
                    p.M(round(xs, prec), round(ys, prec))
                    p.L(round(xn, prec), round(yn, prec))
                # Draw nametag
                tx, ty = (xs+xn)/2, (ys+yn)/2 # point in the middle on straight line
                dt = dom_txt(dname, tx, ty, ux, uy)
                if dname:
                    cx, cy = xs + ux*scale/2, ys + uy*scale/2
                    while len(tseq):
                        circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                        circles.append(nuc_txt(tseq.pop(0), cx, cy, ux, uy))
                        cx, cy = cx + ux*scale, cy + uy*scale
                xs, ys = xn, yn
                yield p, dt
                yield circles

        else: # not enough space, let's draw a rectangle!
            assert slength > linelen
            #detour = get_detour((x0, y0), (x1, y1), slength)
            detour = get_detour_circ((x0, y0), (x1, y1), slength)
            circles = []
            (xs, ys), di = detour[0], 1
            off = scale/2
            paths = []
            for (dname, dlen, dcolor, seq, end) in zip(ns, ls, cs, ss, es):
                tseq = list(seq)
                rsdlen = scale * dlen # remaining scaled domain length
                tlabel, tl = scale * dlen / 2, True
                p = dom_path(dcolor)
                p.M(round(xs, prec), round(ys, prec))
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
                    if a_gt_b(off, seglen): # a >= b
                        off = off - seglen
                        cx, cy = None, None
                    else:
                        cx, cy = xs + ux*off, ys + uy*off
                    # Finish the path if segment length >= remaining domain length
                    if a_eq_b(rsdlen, seglen) or seglen > rsdlen:
                        xn, yn = xs + ux*rsdlen, ys + uy*rsdlen
                        # Draw domain
                        #if end == 'p5':
                        #    p = draw_5prime(p, xs, ys, ux, uy)
                        #    p.L(round(xn, prec), round(yn, prec))
                        #    end = False
                        #elif end == 'p3':
                        if end == 'p3':
                            p.M(round(xs, prec), round(ys, prec))
                            p = draw_3prime(p, xn, yn, ux, uy)
                            end = False
                        else:
                            #p.L(round(xs, prec), round(ys, prec))
                            p.L(round(xn, prec), round(yn, prec))
                        # Draw nucleotides
                        #for nuc in range(dlen):
                        while len(tseq):
                            assert cx is not None
                            circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                            circles.append(nuc_txt(tseq.pop(0), cx, cy, ux, uy),)
                            cx, cy = cx + ux*scale, cy + uy*scale
                            if a_gt_b(distance((xs, ys), (cx, cy)), rsdlen): # a >= b
                                off = scale/2
                                break
                        if a_eq_b(rsdlen, seglen):
                            di+=1
                        xs, ys = xn, yn
                        break
                    # Draw domain
                    #if end == 'p5':
                    #    p = draw_5prime(p, xs, ys, ux, uy)
                    #    p.L(round(xn, prec), round(yn, prec))
                    #    end = False
                    #else:
                        #p.L(round(xs, prec), round(ys, prec))
                    p.L(round(xn, prec), round(yn, prec))
                    # Draw nucleotides
                    while len(tseq):
                        if cx is None:
                            break
                        circles.append(nuc_circle(cx, cy, dcolor, ux, uy))
                        circles.append(nuc_txt(tseq.pop(0), cx, cy, ux, uy))
                        cx, cy = cx + ux*scale, cy + uy*scale
                        if a_gt_b(distance((xs, ys), (cx, cy)), seglen): # a >= b
                            off = distance((xs, ys), (cx, cy)) - seglen
                            break
                    rsdlen -= seglen
                    if tl:
                        tlabel -= seglen
                    di += 1 
                    xs, ys = xn, yn

                yield dt,
                paths.append(p)

            # very dirty
            while di < len(detour):
                xn, yn = detour[di]
                p.L(round(xn, prec), round(yn, prec))
                di += 1
            yield *paths, *circles

