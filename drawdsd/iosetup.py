
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

def draw_complex(stable, ptable, name = '', **kwargs):
    svgC = get_svg_components(stable, ptable, **kwargs)
    return get_drawing(svgC, name)

def main():
    # An minimal workflow using the dsdobjects library.
    import sys

    if False:
        # Example input (pil format):
        pil = """
        sequence x = GGGGGNNNNN: 10
        sequence b = ACGTTNNNNN: 10
        sequence t = 5

        F1 = x( t( b + ) ) t*
        F2 = x( t( b ) ) t*
        """
    else:
        pil = "".join(sys.stdin.readlines())

    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.

    info = read_pil(pil)
    domains = info['domains']
    complexes = info['complexes']

    set_domain_colors(domains)

    for n, c in complexes.items():
        print(f'Drawing complex_{n}')
        stable = list(c.strand_table)
        ptable = list(c.pair_table)
        svg = draw_complex(stable, ptable, name = n)
        svg.save_png(f'complex_{n}.png')


if __name__ == '__main__':
    main()
