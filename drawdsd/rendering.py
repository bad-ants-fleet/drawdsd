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
font_weight = 'bold' # TODO: need some interfacd for plotting parameters.
font_family = 'Helvetica'
domain_stroke_width = 9
backbone_stroke_width = 16
backbone_color = 'gray'
text_font_size = 14
circular_loops = True

# Layout
def get_rgb_palette(num):
    """The default color scheme. """
    num = int(360/(num+1))
    palette = [hls_to_rgb(angle/360, .7, .8) for angle in range(0, 360, num)]
    return [(int(r*255), int(g*255), int(b*255)) for (r,g,b) in palette]

def dom_path(dcolor, sw = None, da = None):
    if sw is None:
        sw = domain_stroke_width
    p = draw.Path(stroke = dcolor, 
                     fill_opacity = 0, 
                     stroke_width = sw, 
                     stroke_dasharray = da,
                     stroke_linecap = 'butt',
                     stroke_linejoin = "round")
    p.dddinfo = 'dom_path'
    return p

def dom_txt(dname, x, y, ux, uy):
    x = x + uy*scale*1.1
    y = y - ux*scale*1.1
    t = draw.Text(dname, text_font_size, 
                     x = round(x, prec), 
                     y = round(y, prec), 
                     font_family = font_family,
                     font_weight = font_weight, 
                     dominant_baseline='middle',
                     text_anchor='middle')
    t.dddinfo = 'dom_txt'
    return t

def pair_shade(x1, y1, x2, y2, x3, y3, x4, y4):
    l = draw.Lines(x1, y1, x2, y2, x3, y3, x4, y4,
                   stroke='none', 
                   fill='lightgray',
                   close='true')
    l.dddinfo = 'pair_shade'
    return l

def h_bonds(x, y, ux, uy):
    p = draw.Path(stroke = 'black',
                  fill_opacity = 0, 
                  stroke_width = 2)
    x0 = x - uy * scale/2
    y0 = y + ux * scale/2
    xn = x - uy * scale*1.2
    yn = y + ux * scale*1.2
    p.M(round(x0, prec), round(y0, prec)).L(round(xn, prec), round(yn, prec))
    p.dddinfo = 'h_bonds'
    return p

def nuc_circle(x, y, c, ux, uy, off = scale/8):
    x = x - uy*off
    y = y + ux*off
    c = draw.Circle(round(x, prec), round(y, prec), round(scale/3, prec), 
                       fill = c, stroke_width = 0, stroke = 'black')
    c.dddinfo = 'nuc_circle'
    return c

def nuc_txt(char, x, y, ux, uy, off = scale/8):
    x = x - uy*off
    y = y + ux*off
    t = draw.Text(char, round(scale*1/2), 
                  round(x, prec), 
                  round(y, prec), 
                  font_family = font_family,
                  font_weight = font_weight, 
                  dominant_baseline='middle',
                  text_anchor='middle')
    t.dddinfo = 'nuc_txt'
    return t

def draw_5prime(p, xs, ys, ux, uy):
    # Look at old commits for different 5' versions ...
    # ... but they all didn't look particularly nice.
    return

def draw_3prime(p, xn, yn, ux, uy):
    # NOTE: The probably most challenging setup for 3' ends 
    # is a 1-nt domain at the end. Improvements are welcome.
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
    if dist == 0:
        raise ValueError(f'Getting distance {dist} for points {p1} and {p2}.')
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

# Main functions.
def draw_stem(i1, i2, i3, i4, p15, p23, p35, p43, 
              color, angle, length, seqT, seqB, nameT, nameB):
    """ Draw stem from initial coordinates. """
    circles = []
    tseqT = list(seqT) if seqT is not None else []
    tseqB = list(seqB) if seqB is not None else []
    # Draw background.
    ps = pair_shade(*i1, *i2, *i3, *i4)
    # Draw top paired domain.
    ux, uy = unitvec(i1, i2)
    dT = dom_path(color)
    dT.M(round(i1[0], prec), round(i1[1], prec))
    if p23:
        dT = draw_3prime(dT, *i2, ux, uy)
    else:
        dT.L(round(i2[0], prec), round(i2[1], prec))
    # Top text layers.
    tx, ty = (i1[0] + i2[0])/2, (i1[1]+i2[1])/2 
    tt = dom_txt(nameT, tx, ty, ux, uy)
    cx, cy = i1[0] + ux*scale/2, i1[1] + uy*scale/2
    while len(tseqT):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt(tseqT.pop(0), cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
    # Draw bottom paired domain.
    ux, uy = unitvec(i3, i4)
    dB = dom_path(color)
    dB.M(round(i3[0], prec), round(i3[1], prec))
    if p43:
        dB = draw_3prime(dB, *i4, ux, uy)
    else:
        dB.L(round(i4[0], prec), round(i4[1], prec))
    # Bottom text layers.
    tx, ty = (i3[0]+i4[0])/2, (i3[1]+i4[1])/2 
    tb = dom_txt(nameB, tx, ty, ux, uy)
    cx, cy = i3[0] + ux*scale/2, i3[1] + uy*scale/2
    while len(tseqB):
        circles.append(h_bonds(cx, cy, ux, uy))
        circles.append(nuc_circle(cx, cy, color, ux, uy))
        circles.append(nuc_txt(tseqB.pop(0), cx, cy, ux, uy))
        cx, cy = cx + ux*scale, cy + uy*scale
    return ps, dT, dB, *circles, tt, tb

def get_circle_data(a, c):
    """Compute the radius and some angles from chord length and arc length.

    Uses scipy.optimize.fsolve for a numerical solution.
    """
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

def get_detour_circ(p0, pN, dtarget):
    """ Return coordinates along a circular backbone.
    """
    dist_0N = distance(p0, pN) # chordlength
    pts = [p0]

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
    points = round(a-0.5) # we don't unscale here.
    assert not a_eq_b(a, c)

    r, phi, psi, off = get_circle_data(a, c)
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
    # shift all by error/2 to center the offset
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

def get_detour_square(p0, pN, dtarget):
    dist_0N = distance(p0, pN)
    ddetour = dtarget - dist_0N
    ux, uy = unitvec(p0, pN, dist_0N)
    rx, ry = vec_rot(ux, uy, -90)
    x1, y1 = p0
    x2, y2 = x1 + rx*ddetour/2, y1 + ry*ddetour/2
    x3, y3 = x2 + ux*dist_0N, y1 + uy*dist_0N
    x4, y4 = x1 + ux*dist_0N, y1 + uy*dist_0N
    assert pN == (x4, y4)
    return [p0, (x2, y2), (x3, y3), pN]
   
def draw_tentacles(dns, dls, dcs, dzs, dos, dss, dms):
    """ Draw loop domains from initial coordinates. 

    TODO: This function would probably benefit from some cleanup.
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
                tseq = list(seq) if seq is not None else []
                sw, da = (None, None) if dname else (2, 5) # fake domain
                p = dom_path(dcolor, sw, da)
                sdlen = scale * dlen
                xn, yn = xs + ux*sdlen, ys + uy*sdlen
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
            if circular_loops:
                detour = get_detour_circ((x0, y0), (x1, y1), slength)
            else:
                detour = get_detour_square((x0, y0), (x1, y1), slength)
            circles = []
            (xs, ys), di = detour[0], 1
            off = scale/2
            paths = []
            for (dname, dlen, dcolor, seq, end) in zip(ns, ls, cs, ss, es):
                tseq = list(seq) if seq is not None else []
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
                        if end == 'p3':
                            p.M(round(xs, prec), round(ys, prec))
                            p = draw_3prime(p, xn, yn, ux, uy)
                            end = False
                        else:
                            p.L(round(xn, prec), round(yn, prec))
                        # Draw nucleotides
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

def estimate_dimensions(svgC):
    """Get the dimensions of the picture from all svg containers. """
    minx, maxx, miny, maxy = 0, 0, 0, 0
    for obj in svgC:
        if obj.dddinfo == 'pair_shade':
            pass
        elif obj.dddinfo == 'nuc_circle':
            pass
        elif obj.dddinfo == 'h_bonds':
            pass
        elif obj.dddinfo == 'nuc_txt':
            pass
        elif obj.dddinfo == 'dom_txt':
            if obj.args['x'] > maxx:
                maxx =  obj.args['x']
            if obj.args['x'] < minx:
                minx =  obj.args['x']
            if obj.args['y'] > maxy:
                maxy =  obj.args['y']
            if obj.args['y'] < miny:
                miny =  obj.args['y']
        elif obj.dddinfo == 'dom_path':
            for x, y in (x[1:].split(',') for x in obj.args['d'].split()):
                if int(float(x)) > maxx:
                    maxx = int(float(x))
                if int(float(x)) < minx:
                    minx = int(float(x))
                if int(float(y)) > maxy:
                    maxy = int(float(y))
                if int(float(y)) < miny:
                    miny = int(float(y))
        else:
            print(f'Unknown object: {x.dddinfo = } {x.TAG_NAME = }.')
    dimx = maxx - minx
    dimy = maxy - miny
    return dimx, dimy, minx, miny

def truncated_domains(segments, truncate = scale/6):
    """Truncate the 5p and/or 3p end of a domain. """
    svg = []
    for segm in segments:
        p = draw.Path(**segm.args)
        x = [[c[0], *c[1:].split(',')] for c in p.args['d'].split()]
        if True: # 5prime
            x1 = float(x[0][1])
            y1 = float(x[0][2])
            x2 = float(x[1][1])
            y2 = float(x[1][2])
            try:
                ux, uy = unitvec((x1, y1), (x2, y2))
                myd = min(distance((x1,y1),(x2,y2)), truncate)
                ps = end_point((x1, y1), (ux, uy), myd)
                x[0][1], x[0][2] = str(round(ps[0], prec)), str(round(ps[1], prec))
            except ValueError:
                # don't truncate if distance between points is too small
                pass
        if True: # 3prime
            xs = float(x[-2][1])
            ys = float(x[-2][2])
            xn = float(x[-1][1])
            yn = float(x[-1][2])
            try:
                ux, uy = unitvec((xn, yn), (xs, ys))
                myd = min(distance((x1,y1),(x2,y2)), truncate)
                ps = end_point((xn, yn), (ux, uy), myd)
                x[-1][1], x[-1][2] = str(round(ps[0], prec)), str(round(ps[1], prec))
            except ValueError:
                # don't truncate if distance between points is too small
                pass
        y = ' '.join([c[0]+','.join(c[1:]) for c in x])
        svg.append(p)
    return svg

def get_backbones(segments):
    """ Returns the backbone from all segments.

    Note: this implementation is not particularly elegant.
    """
    def get_path_ends(path):
        ms = path.args['d'].split()
        assert ms[0][0] == 'M'
        p0 = tuple(ms[0][1:].split(','))
        assert ms[-1][0] == 'L'
        pN = tuple(ms[-1][1:].split(','))
        return p0, pN

    sp0s = {}
    spNs = {}
    segd = {}
    for s in segments:
        p0, pN = get_path_ends(s)
        segd[(p0, pN)] = s
        sp0s[p0] = sp0s.get(p0, []) + [s]
        spNs[pN] = sp0s.get(pN, []) + [s]

    paths = []
    while len(segments):
        seg = segments.pop(0)
        path = [seg]
        p0, pN = get_path_ends(seg)
        seen = []
        while len(segments):
            nseg = segments.pop(0)
            np0, npN = get_path_ends(nseg)
            if npN == p0:
                path.insert(0, nseg)
                p0 = np0
                seen = []
            elif np0 == pN:
                path.append(nseg)
                pN = npN
                seen = []
            else:
                segments.append(nseg)
                seen.append(nseg)
            if len(seen) == len(segments):
                break
        paths.append(path)

    svg = []
    for path in paths:
        p = draw.Path(**(path[0]).args)
        p.args['stroke'] = backbone_color
        p.args['stroke-width'] = backbone_stroke_width
        p.args['stroke-dasharray'] = None,
        p.args['stroke-linecap'] = 'butt'
        d = p.args['d']
        for seg in path[1:]:
            x = seg.args['d'].split()
            y = ' '.join(x[1:])
            p.args['d'] += ' '+y
        svg.append(p)
    return svg

def get_drawing(svgC, name = ''):
    """Get the final svg image from all containers.
    """
    # Initialize the SVG image & background
    dimx, dimy, minx, miny = estimate_dimensions(svgC)

    # Add a 50px frame around the molecule
    dimx += 100
    dimy += 100
    minx -= 50
    miny -= 50
    
    svg = draw.Drawing(dimx, dimy, origin = (minx, miny), displayInline = False)
    svg.append(draw.Rectangle(minx, miny, dimx, dimy, fill = 'white'))
    svg.append(draw.Rectangle(minx+20, miny+20, dimx-40, dimy-40, 
                              stroke_width = 5, stroke = 'darkgray', fill = 'none'))
    
    # Resolve overlap conflicts by collecting similar 
    # objects into separate layers.
    bg_layer = []
    backbone = []
    nuc_layer = []
    txt_layer = []
    for x in svgC:
        if x.dddinfo == 'dom_path':
            backbone.append(x)
        elif x.dddinfo == 'dom_txt':
            txt_layer.append(x)
        elif x.dddinfo == 'pair_shade':
            bg_layer.append(x)
        elif x.dddinfo == 'h_bonds':
            nuc_layer.append(x)
        elif x.dddinfo == 'nuc_circle':
            nuc_layer.append(x)
        elif x.dddinfo == 'nuc_txt':
            txt_layer.append(x)
        else:
            print(f'Unknown object: {x.dddinfo = } {x.TAG_NAME = }.')


    # We truncate domains consistently at the 3' and 5' ends 
    # to give better contrast.
    dom_layer = truncated_domains(backbone)
    bbn_layer = get_backbones(backbone)

    svg.extend(bg_layer)
    svg.extend(bbn_layer)
    svg.extend(dom_layer)
    svg.extend(nuc_layer)
    svg.extend(txt_layer)

    # Draw name in right upper corner
    if name:
        svg.append(draw.Text(name, text_font_size, 
                             x = minx + 30, y = miny + 40,
                             font_family = font_family,
                             font_weight = font_weight))

    return svg


