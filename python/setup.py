#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for installing the model.

To install, run:

    python setup.py install

"""

# Modified from https://github.com/pypa/sampleproject/blob/master/setup.py

from setuptools import setup, find_packages

from codecs import open
from os import path
import sys


if sys.argv[-1] == 'setup.py':
    print("To install lineagesim, run 'python setup.py install'\n")

if sys.version_info[:2] < (2, 7):
    print("hankshaw requires Python 2.7 or later (%d.%d detected)." % sys.version_info[:2])
    sys.exit(-1)


here = path.abspath(path.dirname(__file__))

setup(
    name='lineagesim',
    version='0.2.1',
    description="Model for the publication 'Competition Between Continuously Evolving Lineages in Asexual Populations'",
    long_description="Simulation model described and used in 'Competition Between Continuously Evolving Lineages in Asexual Populations' by Noah Ribeck, Joseph S. Mulka, Luis Zaman, Brian D. Connelly, and Richard E. Lenski. In each simulation, a population grows, mutates, and evolves by selection.",
    url='https://github.com/kerrlab/lineagesim',
    author='Brian D. Connelly',
    author_email='bdc@bconnelly.net',
    license='BSD',

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Artificial Life',
        'License :: OSI Approved :: BSD License'
    ],

    keywords='evolution simulation population mutation interference',

    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    install_requires=['six>=1.10.0', 'numpy>=1.11.0', 'scipy>=0.17.0', 'python-igraph>=0.7.1'],

    extras_require={},
    package_data={},
    data_files=[],

    include_package_data=True,

    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'lineagesim=run_lineagesim:run_simulation',
        ],
    },
)

