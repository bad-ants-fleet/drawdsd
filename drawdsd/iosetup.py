#!/usr/bin/env python

import sys, argparse
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line

from .drawdsd import get_svg_components
from .rendering import get_rgb_palette, get_drawing

def set_domain_lengths(pil_like):
    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.
    pd = read_pil(pil_like)
    for n, d in pd['domains'].items():
        if d.sequence is None:
            d.sequence = 'N' * d.length
            (~d).sequence = 'N' * d.length
    return pd['domains']
 
def set_domain_colors(d):
    palette = get_rgb_palette(len(d))
    for name, dom in sorted(d.items()):
        col = palette.pop(0)
        dom.color = f'rgb{col}'
        (~dom).color = f'rgb{col}'
    return

def draw_complex(cplx, **kwargs):
    stable = list(cplx.strand_table)
    ptable = list(cplx.pair_table)
    svgC, pa, ll, la = get_svg_components(stable, ptable, **kwargs)
    return get_drawing(svgC, cplx.name), pa, ll, la

def main():
    """A wrapper to call drawDSD from the command line."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="A wrapper to call drawDSD from the command line.",
    )
    parser.add_argument("-v", "--verbose", 
            action = "count", default = 0, 
            help = "Verbose output")
    parser.add_argument("-a", "--pair-angles", 
            nargs='+', metavar='<int>:<flt>', default = [],
            help="Specify pair angles, e.g.: -a 1:90 3:180")
    parser.add_argument("-l", "--loop-lengths", 
            nargs='+', metavar='<int>,<int>:<flt>', default = [],
            help="Specify loop lengths, e.g.: -l 0,1:5 1,4:3")
    parser.add_argument("-k", "--loop-angles", 
            nargs='+', metavar='<int>,<int>:<flt>', default = [],
            help="Specify loop angles, e.g.: -l 0,1:-45 1,4:0")
    parser.add_argument("-r", "--rotate", default = 0, type=float,
            help = "Rotate the whole complex by angle.")
    parser.add_argument("-s", "--spacing", default = 0, type=float,
            help = "Include s spacing for 0-loops.")
    args = parser.parse_args()

    pad = {}
    for (k, v) in (x.split(':') for x in args.pair_angles):
        pad[int(k)] = float(v)

    lld = {}
    for (k, v) in (x.split(':') for x in args.loop_lengths):
        sid, did = k.split(',')
        lld[(int(sid), int(did))] = float(v)

    lad = {}
    for (k, v) in (x.split(':') for x in args.loop_angles):
        sid, did = k.split(',')
        lad[(int(sid), int(did))] = float(v)

    # An minimal workflow using the dsdobjects library.
    pil = ''
    if sys.stdout.isatty():
        print("# Reading pil format from STDIN, Use Ctrl-D to finish.")
        pil = "".join(sys.stdin.readlines())
    if not pil:
        pil = """
        sequence a = NNNNNNNNNN: 10
        sequence t = NNNNN: 5
        sequence x = NNNNNNNNNN: 10
        sequence b = NNNNNNNNNN: 10
        F1 = x( t( b + ) ) t*
        """
        print(f"# Using example pil input:\n{pil}")

    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.

    info = read_pil(pil)
    domains = info['domains']
    complexes = info['complexes']
    set_domain_colors(domains)

    for n, cplx in complexes.items():
        print(f'Drawing complex {n}:')
        if not (pad or lld or lad):
            svg, pa, ll, la = draw_complex(cplx, rotate = args.rotate, spacing = args.spacing)
            print(f' Plotting parameters:')
            print(f'  - pair-angles = {pa}')
            print(f'  - loop-lengths = {ll}')
            print(f'  - loop-angles = {la}')
            svg.save_png(f'complex_{n}.png')
        else:
            _, pa, ll, la = draw_complex(cplx)
            for k, v in pad.items():
                pa[k] = v
            for (k, l), v in lld.items():
                ll[k][l] = v 
            for (k, l), v in lad.items():
                la[k][l] = v 
            svg, pa, ll, la = draw_complex(cplx, 
                                       pair_angles = pa, 
                                       loop_lengths = ll, 
                                       loop_angles = la, 
                                       rotate = args.rotate, 
                                       spacing = args.spacing)
            print(f' Plotting parameters:')
            print(f'  - pair-angles = {pa}')
            print(f'  - loop-lengths = {ll}')
            print(f'  - loop-angles = {la}')
            svg.save_png(f'complex_{n}.png')

if __name__ == '__main__':
    main()
