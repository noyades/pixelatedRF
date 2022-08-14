#!/usr/bin/env python
from setuptools import setup

"""
with open("README.rst") as fin:
    long_description = fin.read()
"""

with open("vt_rrfc/__init__.py") as fin:
  for line in fin:
    if line.startswith("__version__ ="):
      version = eval(line[14:])
      break

setup(name='vt_rrfc',
      version=version,
      description='Tools for generating gds for inverse design in PCB and IC',
      author='Jeff Walling',
      author_email='jswalling@vt.edu',
      install_requires=['gdspy','numpy','scipy'],
      url="https://github.com/noyades/pixelatedRF",
      packages=['vt_rrfc']
      )
