# DrawDSD -- SVG images of domain-level complexes 

Tired of drawing domain-level complexes by hand? Tired of automatically
generated images that are *almost* exactly like you want them? Then try this.

## Test and install
```sh
    python -m pytest tests/ -vs
    python setup.py install
```

## Quickstart

### 1. Using Jupyter notebooks:
**This interface is in prepration, it does not work yet!** 

### 2. Customizing the SVG output:
This example assumes you are somewhat comfortable working with the [dsdobjects]
library (although you may not need more than the few lines of code below).
Customization involves three types of parameters:
 - the color of domains
 - the angles of paired domains
 - the distances of unpaired domains

```py
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line
from drawdsd import draw_complex, get_default_plot_params
from drawdsd.rendering import get_drawing, get_rgb_palette

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
mycplx = read_pil_line('A = r b( g r b( l ) y r b( g + l ) y r b( g + l ) y l ) y')

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
svg.savePng(f'drawing_complex_{mycplx.name}.png')
```
You can find a modified version of this script in the [examples] folder: [customize.py].

## Version
 - v0.1

## TODOs (ranking from high to low priority)
~~- draw 5' squares and 3' arrow-heads (test for 'A = x( a ) y( a )')~~
 - force text on top of drawings.
 - support "p3 attachment" (e.g. for 'A = x( r + r( + ) x*( + ) t2*( + ) )')
 - provide interface for pair-lengths and loop-angles.
 - provide interface for plotting parameters such as path width, font size, etc.
 - rotate complexes and fix alignment of domain names (e.g. all pair_angles +5)
 - find nice examples for rendering and make a notebook gallery.
 - draw sequence information.
 - provide better default rendering (e.g naview-like, forgi-like)
 - provide additional high-level input from PIL file only.
 - support specification of pair-angles and loop-lengths in kernel string.
 - support Bezier curves for loops (as alternative to rectangles).

# Have fun!

[dsdobjects]: <https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects>
[peppercornenumerator]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[examples]: <https://github.com/bad-ants-fleet/drawdsd/tree/master/examples>
[customize.py]: <https://github.com/bad-ants-fleet/drawdsd/tree/master/examples/customize.py>
