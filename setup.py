#!/usr/bin/env python
import os, sys
#from numpy.distutils.core import setup, Extension
from setuptools import find_packages
from distutils.core import setup

if sys.version_info.major != 3:
    sys.exit("scikit-ribo can only be used with Python 3. You are currently "
             "running Python %d." % sys.version_info.major)

setup(
    name='scikit_ribo',
    version='0.2.4b1',
    description = 'A scikit framework for joint analysis of Riboseq and RNAseq data',
    long_description=open('pypi_readme.rst').read(),
    url="https://github.com/hanfang/scikit-ribo",
    author = 'Han Fang',
    author_email = 'hanfang.cshl@gmail.com',
    license='GPLv2',
    scripts=['scikit_ribo/scikit-ribo-run.py','scikit_ribo/scikit-ribo-build.py','scikit_ribo/plot_ribo.py'],
    packages=find_packages(),
    install_requires=[
                      'colorama>=0.3.7',
                      'gffutils>=0.8.7.1',
                      'glmnet_py>=0.1.0b2',
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
        'Operating System :: Unix',
                ],
    keywords='bioinformatics genomics glm glmnet ridge riboseq',
        )
