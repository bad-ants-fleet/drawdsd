#
# Produce an SVG file from data.
#
import logging
log = logging.getLogger(__name__)

import numpy as np
import drawsvg as draw
from colorsys import hls_to_rgb

dfont = 'bold' # TODO: need some interfacd for plotting parameters.
ffamily = 'Courier'

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
    dimx = maxx - minx
    dimy = maxy - miny
    return dimx, dimy, minx, miny

def get_rgb_palette(num):
    num = int(360/(num+1))
    palette = [hls_to_rgb(angle/360, .6, .8) for angle in range(0, 360, num)]
    return [(int(r*255), int(g*255), int(b*255)) for (r,g,b) in palette]

def draw_stem(i1, i2, i3, i4, p15, p23, p35, p43, color, angle, length, nameT, nameB):
    bground = draw.Path(stroke_width=0, fill='black', fill_opacity=.3)
    dompath = draw.Path(stroke_width=5, stroke=color, fill_opacity=0)

    (x1, y1) = i1
    (x10, y10) = i1
    (x11, y11) = i1
    (x12, y12) = i1
    (x2, y2) = i2
    (x3, y3) = i3
    (x30, y30) = i3
    (x31, y31) = i3
    (x32, y32) = i3
    (x4, y4) = i4

    bground.M(round(x1), round(y1)).L(
              round(x2), round(y2)).L(
              round(x3), round(y3)).L(
              round(x4), round(y4))

    oS = 6
    if p15:
        x10, y10 = x1+(x2-x1)/length*(oS+5), y1+(y2-y1)/length*(oS+5)
        x11, y11 = x10-(y2-y1)/length*oS, y10+(x2-x1)/length*oS
        x12, y12 = x11-(x2-x1)/length*oS, y11-(y2-y1)/length*oS
        x1,  y1  = x12+(y2-y1)/length*oS, y12-(x2-x1)/length*oS, 
    if p35:
        x30, y30 = x3+(x4-x3)/length*(oS+5), y3+(y4-y3)/length*(oS+5)
        x31, y31 = x30-(y4-y3)/length*oS, y30+(x4-x3)/length*oS
        x32, y32 = x31-(x4-x3)/length*oS, y31-(y4-y3)/length*oS
        x3, y3 = x32+(y4-y3)/length*oS, y32-(x4-x3)/length*oS, 

    dompath.M(round(x10), round(y10)).L(
              round(x11), round(y11)).L(
              round(x12), round(y12)).L(
              round(x1), round(y1)).L(
              round(x2), round(y2)).M(
              round(x30), round(y30)).L(
              round(x31), round(y31)).L(
              round(x32), round(y32)).L(
              round(x3), round(y3)).L(
              round(x4), round(y4))

    cT = 10
    tx, ty = (x1+x2)/2, (y1+y2)/2 # point in the middle on straight line
    tx, ty = tx-(y2-y1)/length*cT, ty+(x2-x1)/length*cT
    tagT = draw.Text(nameT, 14, x = round(tx), y = round(ty), 
                        font_family = ffamily,
                        font_weight = dfont, 
                        dominant_baseline='middle',
                        text_anchor='middle')

    tx, ty = (x3+x4)/2, (y3+y4)/2 # point in the middle on straight line
    tx, ty = tx-(y4-y3)/length*cT, ty+(x4-x3)/length*cT
    tagB = draw.Text(nameB, 14, x = round(tx), y = round(ty), 
                        font_family = ffamily,
                        font_weight = dfont, 
                        dominant_baseline='middle',
                        text_anchor='middle')

    return bground, dompath, tagT, tagB

def draw_tentacles(dns, dls, dcs, dzs, dos, dms):
    # names, lengths, colors, startpoint, endpoint
    for ns, ls, cs, (x0, y0), (x1, y1), es in zip(dns, dls, dcs, dzs, dos, dms):
        length = sum(ls) # total length of concatenated domains
        linelen = np.sqrt((x1-x0)**2 + (y1-y0)**2) # end-to-end distance

        # Append a fake domain if the distance between points is 
        # longer than the length of all domains.
        if not np.isclose(length, linelen) and linelen > length:
            ns.append('')
            ls.append(linelen-length)
            cs.append('black')
            length = sum(ls)

        if np.isclose(length, linelen):
            xs, ys = x0, y0 # start coordinates
            for (dname, dlen, dcolor, end) in zip(ns, ls, cs, es):
                # Draw path
                (sw, da) = (5, None) if dname else (2, 5) # fake domain
                p = draw.Path(stroke = dcolor, fill_opacity = 0,
                              stroke_width = sw, stroke_dasharray = da)
                xn, yn = xs+(x1-x0)/length*dlen, ys+(y1-y0)/length*dlen

                if end == 'p5':
                    oS = 6
                    ux = (x1-x0)/length # unit x component
                    uy = (y1-y0)/length # unit y component
                    x50, y50 = xs+ux*(oS+5), ys+uy*(oS+5)
                    x51, y51 = x50-uy*oS, y50+ux*oS
                    x52, y52 = x51-ux*oS, y51-uy*oS
                    xs, ys = x52+uy*oS, y52-ux*oS
                    p.m(round(x50), round(y50)).L(
                        round(x51), round(y51)).L(
                        round(x52), round(y52)).L(
                        round(xs), round(ys)).L(
                        round(xn), round(yn))#.L( round(x3),round(y3))
                elif end == 'p3':
                    oS = 8
                    ux = (x1-x0)/length # unit x component
                    uy = (y1-y0)/length # unit y component
                    x30, y30 = xn-ux*(oS+5), yn-uy*(oS+5)
                    x32, y32 = x30-uy*oS, y30+ux*oS
                    p.m(round(xs), round(ys)).L(
                        round(x30), round(y30)).L(
                        round(xn), round(yn)).L(
                        round(x32),round(y32))
                else:
                    p.m(round(xs), round(ys)).L(round(xn), round(yn))

                # Draw nametag
                cT = 10
                tx, ty = (xs+xn)/2, (ys+yn)/2 # point in the middle on straight line
                tx, ty = tx-(yn-ys)/dlen*cT, ty+(xn-xs)/dlen*cT
                t = draw.Text(dname, 14, x = round(tx), y = round(ty), 
                        font_family = ffamily,
                        font_weight = dfont, 
                        dominant_baseline='middle',
                        text_anchor='middle')
                #t = draw.Circle(round(tx), round(ty), 10, 
                #            fill='white', stroke_width=2, stroke='black')
                yield p, t
                xs, ys = xn, yn

        else: # not enough space, let's draw a rectangle!
            assert length > linelen
            missing = length - linelen
            # The points of the rectangle: M(x0, y0).L(rx2, ry2).L(rx3, ry3).L(x1, y1)
            (rx2, ry2) = x0-(y1-y0)/linelen*missing/2, y0+(x1-x0)/linelen*missing/2
            (rx3, ry3) = x1-(y1-y0)/linelen*missing/2, y1+(x1-x0)/linelen*missing/2
            slen = np.sqrt((rx2-x0)**2 + (ry2-y0)**2)
            assert np.isclose(length, linelen+2*slen)

            xs, ys = x0, y0
            cseg, ci = [(rx2, ry2), (rx3, ry3), (x1, y1)], 0
            for (dname, dlen, dcolor) in zip(ns, ls, cs):
                tl = dlen/2
                p = draw.Path(stroke = dcolor, fill_opacity = 0, stroke_width = 5,
                        stroke_linejoin="round")
                p.M(round(xs), round(ys)) # move to start.
                while True:
                    xn, yn = cseg[ci]
                    plen = np.sqrt((xn-xs)**2 + (yn-ys)**2) # Rename

                    # Draw nametag if we are further than half-way.
                    if np.isclose(tl, plen) or tl < plen:
                        tx, ty = xs+(xn-xs)/plen*tl, ys+(yn-ys)/plen*tl

                        cT = 20 if xs > xn else 10
                        tx, ty = tx-(yn-ys)/plen*cT, ty+(xn-xs)/plen*cT
                        t = draw.Text(dname, 14, x = round(tx), y = round(ty), font_weight = dfont, 
                                      text_anchor='middle', valign='center')
                        tl = 2 * dlen
                    else:
                        tl -= plen

                    # Finish the path if pathlength >= domainlength
                    if np.isclose(dlen, plen) or plen > dlen:
                        xn, yn = xs+(xn-xs)/plen*dlen, ys+(yn-ys)/plen*dlen
                        p.L(round(xn), round(yn))
                        xs, ys = xn, yn
                        break
                    p.L(round(xn), round(yn))
                    dlen -= plen
                    ci += 1 
                    xs, ys = xn, yn
                yield p, t

