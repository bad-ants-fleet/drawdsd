# DrawDSD -- SVG images of domain-level complexes 

Tired of drawing domain-level complexes by hand? Tired of automatically
generated images that are *almost* exactly like you want them? Then try this.

## Test and install
```sh
    python -m pytest tests/ -vs
    python setup.py install
```

## Quickstart

### 1. Using jupyter notebooks:
**This interface is in prepration, it does not work yet!** 
Instead, call the preliminary main function from the command line ...
```
    python drawdsd/drawdsd.py
```
... and then look at the file *example.png*.


### 2. A high-level interface:
**This interface is in prepration, it does not work yet!** 
This interface does not support customization, so it might not be what you are
looking for, but it is a good place to start. The input is provided in form
of a pilfile (e.g. see the file format used for [peppercornenumerator]).

```py
import drawSvg as draw
from drawdsd import draw_kernel

# Let's get started with a white backround.
d = draw.Drawing(800, 800, origin = (-250, -550), displayInline = False)
d.append(draw.Rectangle(-100, -400, 500, 500, fill='white'))

# Specify a complex in kernel format.
pilstring = """
length a = 7
length b = 5
length c = 15
length x = 5
length y = 15
length z = 15

A = a x( b c y( a* + c ) z( b* ) )
"""

# Draw the complex on the canvas.
for (svgs, pa, ll) in draw_kernel(pilstring, is_file = False):
    d.extend(svgs)

# Well, let's look it!
d.savePng('example.png')
```

### 3. Customizing the SVG output:
This example assumes you are somewhat comfortable working with the [dsdobjects]
library (altough you may not need more than the few lines of code below).
Customization involves three types of parameters:
 - the color of domains
 - the angles of paired domains
 - the distances of unpaired domains

```py
from dsdobjects.objectio import set_io_objects, read_pil, read_pil_line
from drawdsd import draw_complex, get_rgb_palette, get_default_plot_params
from drawdsd.rendering import get_drawing

set_io_objects() # Using the default Domain, Complex objects of dsdobjects.

# Let's separate the initialization of Domain objects ...
read_pil('''
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
ll[0][1] = 20 # Force the distance of loop #1 (2nd) on strand #0 (1st).
ll[1][1] = 20

# Third, get the SVG objects of the complex!
svgC, pa, ll = draw_complex(mycplx, pair_angles = pa, loop_lengths = ll)

# Last, draw the complex!
svg = get_drawing(svgC)
svg.append(draw.Text(f'{mycplx.name}:', 14, x = -25, y = 0, 
                     font_weight = 'bold', text_anchor='middle', valign='center'))
svg.savePng(f'drawing_complex_{mycplx.name}.png')
```

## Version
 - v0.1

## TODOs
 - draw 5' and 3' ends.
 - draw sequence information.
 - support Bezier curves for loops (as alternative to rectangles).
 - support specification of pair-angles and loop-lengths in kernel string.
 - provide interface for pair-lengths and loop-angles.
 - provide interfaces as shown in README

# Have fun!

[dsdobjects]: <https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects>
[peppercornenumerator]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
