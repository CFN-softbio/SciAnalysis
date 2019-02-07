#!/usr/bin/env python
import numpy as np

# from distutils.core import setup
import setuptools

setuptools.setup(name='SciAnalysis',
                 author='Kevin Yager',
                 packages=setuptools.find_packages(),
                 description="SciAnalysis scripts for processing image files.",
                 include_dirs=[np.get_include()],
                 author_email='kyager@bnl.gov',
                 keywords='Analysis',
                 )
