#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    LONG_DESCRIPTION = fh.read()

setup(
    name = 'drawdsd',
    version = '0.1',
    description = 'SVG images of domain-level complexes',
    long_description = LONG_DESCRIPTION,
    long_description_content_type = "text/markdown",
    author = 'Stefan Badelt',
    author_email = 'bad-ants-fleet@posteo.eu',
    maintainer = 'Stefan Badelt',
    maintainer_email = 'bad-ants-fleet@posteo.eu',
    url = 'https://github.com/bad-ants-fleet/drawdsd',
    license = 'MIT',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Programming Language :: Python :: 3.9',
        'Intended Audience :: Science/Research',
        ],
    python_requires = '>=3.9',
    install_requires = ['drawSvg',
                        'dsdobjects'],
    packages = find_packages(),
    test_suite = 'tests',
)

