# drawdsd -- SVG images of domain-level complexes 

Tired of drawing domain-level complexes by hand? Tired of automatically
generated images that are *almost* exactly like you want them? Try this.

## Install
```sh
    pip install git+https://github.com/bad-ants-fleet/drawdsd.git
```

## Quickstart
A good starting point is the [examples] directory. You can find a Jupyter
notebook [examples.ipynb] to learn about the most important drawing features.
There are also input files that can be visualized using the command line tool
*drawDSD*. It is probably worth looking at the code of this script if you want
get started with your own code based on drawdsd. The main function can be found 
in the file [iosetup.py].

**Note**: For convenience, the drawdsd library uses [dsdobjects] to translate 
text-based "kernel notation" input into Python objects. An introduction about 
kernel notation can be found at the [peppercornenumerator] project. Importantly,
the objects provided by [dsdobjects] are immutable with respect to certain 
properties, e.g. the length of a domain cannot be changed, after it has been 
initialized. 

## Version
 - v0.3: major code update, can draw all pk-free structures, new colors, 3'ends, backbone, etc.

### TODO:
 - Provide visualization of unpaired strands. (Currently, an input complex needs at least one paired domain.)
 - support pseudoknots.
 - support simple reaction scheme.

[dsdobjects]: <https://github.com/DNA-and-Natural-Algorithms-Group/dsdobjects>
[peppercornenumerator]: <https://github.com/DNA-and-Natural-Algorithms-Group/peppercornenumerator>
[examples]: <https://github.com/bad-ants-fleet/drawdsd/tree/master/examples>
[examples.ipynb]: <https://github.com/bad-ants-fleet/drawdsd/blob/master/examples/examples.ipynb>
[iosetup.py]: <https://github.com/bad-ants-fleet/drawdsd/blob/master/drawdsd/iosetup.py>
