#!/usr/bin/env python
import os, sys
#from numpy.distutils.core import setup, Extension
from setuptools import find_packages

from distutils.core import setup

setup(
    name='scikit_ribo',
    version='0.2.0b11',
    description = 'A scikit framework for joint analysis of Riboseq and RNAseq data',
    long_description=open('README.md').read(),
    url="http://pypi.python.org/pypi/scikit-ribo/",
    author = 'Han Fang',
    author_email = 'hanfang.cshl@gmail.com',
    license='GPLv2',
    scripts=['scikit_ribo/scikit-ribo-run.py','scikit_ribo/scikit-ribo-build.py'],
    packages=find_packages(),
    install_requires=['conda>=4.2.13',
                      'colorama>=0.3.7',
                      'gffutils>=0.8.7.1',
                      'joblib>=0.10.3',
                      'matplotlib>=1.5.1',
                      'numpy>=1.11.2',
                      'pandas>=0.19.2',
                      'pybedtools>=0.7.8',
                      'pyfiglet>=0.7.5',
                      'pysam>=0.9.1.4',
                      'scikit-learn>=0.18',
                      'scipy>=0.18.1',
                      'seaborn>=0.7.0',
                      'termcolor>=1.1.0',
                      ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
                ],
    keywords='bioinformatics genomics glm glmnet ridge riboseq',
        )
