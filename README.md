# DrawDSD -- SVG images of domain-level complexes 

Tired of drawing domain-level complexes by hand? Tired of automatically
generated images that are *almost* exactly like you want them? Then try this.

## Test and install
```sh
    python -m pytest tests/ -vs
    pip install .
```

## Quickstart

### 1. Using Jupyter notebooks:
see examples/examples.ipynb

### 2. Customizing the SVG output:
Here is a small example script using the default rendering. The
pair angles, loop lengths, and loop angles are adjustable, see
the notebook for details.

```py
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
"""

set_io_objects() # Using the default Domain, Complex objects of dsdobjects.

info = read_pil(pil)
domains = info['domains']
complexes = info['complexes']

set_domain_colors(domains)

for n, cplx in complexes.items():
    print(f'Drawing complex_{n}')
    svg, pa, ll, la = draw_complex(cplx)
    print(f'pair-angles = {pa}')
    print(f'loop-lengths = {ll}')
    print(f'loop-angles = {la}')
    svg.save_png(f'complex_{n}.png')

```

## Version
 - v0.3: 

## TODO:
 - provide visualization of single strands?

[dsdobjects]: <https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects>
[peppercornenumerator]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[examples]: <https://github.com/bad-ants-fleet/drawdsd/tree/master/examples>
[customize.py]: <https://github.com/bad-ants-fleet/drawdsd/tree/master/examples/customize.py>
