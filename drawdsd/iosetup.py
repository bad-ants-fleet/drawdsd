
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line

from .drawdsd import get_default_plot_params, draw_complex
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

#def get_plot_params(kernel, pa = None, ll = None, seq = None):
#    cplx = read_pil_line(kernel)
#    return get_default_plot_params(cplx)
#
#def draw(kernel, pa = None, ll = None, seq = None, direction = 1):
#    cplx = read_pil_line(kernel)
#    svgC = draw_complex(cplx, pair_angles = pa, loop_lengths = ll, sequence = seq)
#    return get_drawing(svgC)


