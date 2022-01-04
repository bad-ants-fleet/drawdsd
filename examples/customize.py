import drawSvg as draw
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line
from drawdsd import draw_complex, get_default_plot_params
from drawdsd.rendering import get_drawing, get_rgb_palette

def main():
    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.
    
    # Let's separate the initialization of Domain objects ...
    _ = read_pil('''
             length a = 7
             length b = 5
             length c = 15
             length x = 5
             length y = 15
             length k = 5
             length r = 5
             length g = 5
             length l = 5
             ''')
    
    # ... from the initialization of a specific complex object via a kernel string:
    #mycplx = read_pil_line('A = k a( b y y c( y + ) ) k')
    mycplx = read_pil_line('A = r b( g r b( l ) y r b( g + l ) y r b( g + l ) y l ) y')
    #mycplx = read_pil_line('A = r b( b( l ) b( g + l ) l ) y')
    
    #
    # Customization 1: Let's choose the color for each domain.
    #
    palette = get_rgb_palette(len(mycplx.domains))
    for dom in mycplx.domains:
        if hasattr(dom, 'color'):
            continue
        col = palette.pop(0)
        dom.color = f'rgb{col}'
        (~dom).color = f'rgb{col}'
            
    #
    # Customization 2 & 3: Let's choose "pair angles" and "loop lengths".
    #
    
    # First, get the defaults.
    pa, ll = get_default_plot_params(mycplx)
    
    # Second, change the defaults. (Values are 0-based and in order of the kernel string input).
    pa[3] = 180   # Force the angle of the paired domain #3 (4th).
    ll[0][1] = 10 # Force the distance of loop #1 (2nd) on strand #0 (1st).
    ll[1][1] = 10
    
    # Third, get the SVG objects of the complex!
    svgC, pa, ll = draw_complex(mycplx, pair_angles = pa, loop_lengths = ll)
    
    # Last, draw the complex!
    svg = get_drawing(svgC)
    svg.append(draw.Text(f'{mycplx.name}:', 14, x = -25, y = 0, 
                         font_weight = 'bold', text_anchor='middle', valign='center'))
    svg.savePng(f'drawing_complex_{mycplx.name}.png')

if __name__ == '__main__':
    main()
