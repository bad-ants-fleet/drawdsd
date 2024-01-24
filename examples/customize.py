#!/usr/bin/env python

import random
import drawsvg as draw
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line
from dsdobjects.iupac_utils import reverse_wc_complement
from drawdsd.iosetup import draw_complex
from drawdsd.drawdsd import get_default_plot_params
from drawdsd.rendering import get_drawing, get_rgb_palette

# A few examples.
kerneldrawings = {
    'x': {'kernel': 'x = b( b*( ) + ) l' ,
            'pa': {},#{0:90, 1: 90, 2:90}, #{0: 90}, #{3: 180},
          'll': {(0,2):22}#{(0,2): 0, (0,3):2 }, # {(0,1):20}
           },
    'y': {'kernel': 'y = x a x L0*( S0 a L0( b( + S1 c* ) ) a*( d*( S2 e ) f +  ) ) a y ',
          'pa': {3:-90, 1:-30},#{0:10, 1:10, 3:-90, 2:40}, #{0: 90}, #{3: 180},
          'll': {(0,1):10}#{(0,1):19, (2,2):None}#{(0,0): 14, (2,2):8}
         },

    'b': {'kernel': 'b = x( r + r( h( k ) a + b ) x*( + ) h*(  ) )',
          'pa': {0:0, 1:0, 2: 0, 3: 0},
          'll': {(0,1):15, (2,1): 0},
          'la': {(0,0): 49, (0,1):90, (1,3): -10, (2,0): -40},
          },

    'A': {'kernel': 'A = x a( y a*( z ) u a( y + x ) u a*( x + ) v ) z',
        #'pa': {},#{0:45, 1:135, 2: 45, 3: -45},#{3: 180},
        'pa': {3: 180},
        'll': {}#{(1,1): 5}
         },
    'A1': {'kernel': 'A1 = x a( y a*( z ) u a( y + x ) u a*( h a*( x + ) ) v ) z',
           'pa': {1: 90, 3: 270, 4:250},
           'll': {(1,2):12.9}
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
          'pa': {1:90},
           'll': {},
         },
    # Feature request: loop angles!
    'D1': {'kernel': 'D1 = a b( b( c ) b( x + c ) l ) y',
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
             length L0 = 15
             length S0 = 15
             length S1 = 15
             length S2 = 5
             length a = 10 
             length b = 10
             length c = 10
             length d = 10
             length e = 1
             length f = 10
             length g = 15
             length h = 15
             length i = 15
             length j = 15
             length k = 15
             length l = 15
             length r = 15
             length u = 5
             length v = 5
             length w = 5
             length x = 15
             length y = 5
             length z = 5
             ''')
    
    # ... from the initialization of a specific complex object via a kernel string:
    mycplx = read_pil_line(ddict['kernel'])
   
    # Customization 1: Let's choose the color for each domain.
    #
    palette = get_rgb_palette(len(mycplx.domains))
    for dom in sorted(mycplx.domains):
        if hasattr(dom, 'color'):
            continue
        col = palette.pop(0)
        dom.color = f'rgb{col}'
        (~dom).color = f'rgb{col}'
        if dom.sequence:
            continue
        dom.sequence = ''.join([random.choice('ACGT') for _ in range(len(dom))])
        (~dom).sequence = reverse_wc_complement(dom.sequence, material = 'DNA')
        #print(dom, dom.sequence)
            
    #
    # Customization 2 & 3: Let's choose "pair angles" and "loop lengths".
    #
    
    stable = list(mycplx.strand_table)
    ptable = list(mycplx.pair_table)

    # First, get the defaults.
    pa, ll, la = get_default_plot_params(stable, ptable)
   
    # Second, change the defaults. All values are 0-based and in order of the
    # kernel string input: 
    for k, v in ddict['pa'].items():
        # Force the angle of the kth paired domain to v.
        pa[k] = v 
    for (k, l), v in ddict['ll'].items():
        # Force the distance of the kth loop on the lth strand to v.
        ll[k][l] = v 

    if 'la' in ddict:
        for (k, l), v in ddict['la'].items():
            # Force the distance of the kth loop on the lth strand to v.
            la[k][l] = v 
   
    # Third, get the SVG objects of the complex!
    svg = draw_complex(mycplx, pair_angles = pa, loop_lengths = ll, loop_angles = la)
    svg.save_png(f'complex_{mycplx.name}.png')
    #svg.save_svg(f'complex_{mycplx.name}.svg')


if __name__ == '__main__':
    main()
