#!/usr/bin/env python

import drawSvg as draw
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line
from drawdsd import draw_complex, get_default_plot_params
from drawdsd.rendering import get_drawing, get_rgb_palette

# A few examples.
kerneldrawings = {
    'A': {'kernel': 'A = x a( y a*( z ) u a( y + x ) u a*( x + ) v ) z',
           'pa': {3: 180},
           'll': {}
         },
    'A1': {'kernel': 'A1 = x a( y a*( z ) u a( y + x ) u a*( x + ) v ) z',
           'pa': {1: 90, 3: 270},
           'll': {}
         },
 
    'B': {'kernel': 'B = z g( b c h( l ) y z b( g + l ) y z b( x + y* ) y z ) y',
           'pa': {0: -45, 1: 45, 2: -45, 3: 180+45},
           'll': {(0, 1): 10}
         },
    # Bug: misplaced domain name!
    'C': {'kernel': 'C = z g( b c h( l ) y z b( g + l ) y z b( g + l ) y l ) y',
           'pa': {0: -45, 1: 45, 2: -45, 3: 180+45},
           'll': {}
         },
    'D': {'kernel': 'D = a b( b( c ) b( x + c ) l ) y',
           'pa': {1: 90},
           'll': {}
         },
    # Feature request: loop angles!
    'D1': {'kernel': 'D = a b( b( c ) b( x + c ) l ) y',
           'pa': {1: 90},
           'll': {(0, 1): 7, (0, 3): 7}
         },
    'E': {'kernel': 'E = x b( c( y + ) ) y* ',
           'pa': {},
           'll': {}
         },
    # Bug (5' and 3' end are connected!)
    'F': {'kernel': 'F = x( a ) y( a )',
        'pa': {0: 180},
           'll': {}
         },
    }

def main():
    """ A playground for trying plots.
    """
    # Choose one of the examples above ...
    ddict = kerneldrawings['A1']

    set_io_objects() # Using the default Domain, Complex objects of dsdobjects.
    
    # Let's separate the initialization of Domain objects ...
    _ = read_pil('''
             length a = 10 
             length b = 10
             length c = 10
             length d = 10
             length e = 10
             length f = 10
             length g = 15
             length h = 15
             length i = 15
             length j = 15
             length k = 15
             length l = 15
             length u = 5
             length v = 5
             length w = 5
             length x = 5
             length y = 5
             length z = 5
             ''')
    
    # ... from the initialization of a specific complex object via a kernel string:
    mycplx = read_pil_line(ddict['kernel'])
   
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
   
    # Second, change the defaults. All values are 0-based and in order of the
    # kernel string input: 
    for k, v in ddict['pa'].items():
        # Force the angle of the kth paired domain to v.
        pa[k] = v 
    for (k, l), v in ddict['ll'].items():
        # Force the distance of the kth loop on the lth strand to v.
        ll[k][l] = v 
   
    # Third, get the SVG objects of the complex!
    svgC, pa, ll = draw_complex(mycplx, pair_angles = pa, loop_lengths = ll)
    
    # Last, draw the complex!
    svg = get_drawing(svgC)
    #svg.append(draw.Text(f'{mycplx.name}:', 14, x = -25, y = 0, 
    #                     font_weight = 'bold', text_anchor='middle', valign='center'))
    #svg.savePng(f'complex_{mycplx.name}.png')
    #svg.saveSvg(f'complex_{mycplx.name}.svg')
    svg.savePng(f'current_complex.png')

if __name__ == '__main__':
    main()
