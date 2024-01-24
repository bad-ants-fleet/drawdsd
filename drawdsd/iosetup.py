#!/usr/bin/env python

import sys
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
    svgC = get_svg_components(stable, ptable, **kwargs)
    return get_drawing(svgC, cplx.name)

def main():
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

        # DrawDSD does not support single strands atm, 
        # so here the toehold is assumed to be paired!
        A = a t( x + )
        B = x t( b + )
        F1 = x( t( b + ) ) t*
        F2 = a t( x( + t* ) )
        F3 = a t( x( t* ) )
        F4 = a t( t x( a ) )
        """
        print(f"# Using example pil input:\n{pil}")

    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.

    info = read_pil(pil)
    domains = info['domains']
    complexes = info['complexes']

    set_domain_colors(domains)

    for n, cplx in complexes.items():
        print(f'Drawing complex_{n}')
        svg = draw_complex(cplx)
        svg.save_png(f'complex_{n}.png')


if __name__ == '__main__':
    main()
