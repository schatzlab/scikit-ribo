![logo](docs/images/logo.png)

## *scikit-ribo* 
[![PyPI version](https://badge.fury.io/py/scikit-ribo.svg)](https://badge.fury.io/py/scikit-ribo)
[![GPL Licence](https://badges.frapsoft.com/os/gpl/gpl.svg?v=103)](https://opensource.org/licenses/GPL-3.0/)
[![Documentation Status](https://readthedocs.org/projects/scikit-ribo/badge/?version=latest)](http://scikit-ribo.readthedocs.io/en/latest/?badge=latest)

#### - Accurate inference and robust modelling of translation dynamics at codon resolution with Riboseq data
#### https://github.com/hanfang/scikit-ribo

--------
## Documentation
Read the Docs: [![Documentation Status](https://readthedocs.org/projects/scikit-ribo/badge/?version=latest)](http://scikit-ribo.readthedocs.io/en/latest/?badge=latest) or click [me](http://scikit-ribo.readthedocs.io/en/latest/)

## Contact

#### Han Fang
#### Stony Brook University & Cold Spring Harbor Laboratory
#### Email: hanfang.cshl@gmail.com

## Requirement: 
#### Environment: Python3, Linux

Recommend setting up your environment with [Conda](https://conda.io/docs/intro.html)

#### Dependencies:

| Dependencies | Version >= |
| ------------- |:-------------:|
| bedtools | 2.26.0 |

When using `pip install scikit-ribo`, all the following dependencies will be pulled and installed automatically.

| Python package| Version >= |
| ------------- |:-------------:|
| colorama | 0.3.7 |
| glmnet-py | 0.1.0b |
| gffutils | 0.8.7.1 |
| matplotlib | 1.5.1 |
| numpy | 1.11.2 |
| pandas | 0.19.2 |
| pybedtools | 0.7.8 | 
| pyfiglet | 0.7.5 | 
| pysam | 0.9.1.4 |
| scikit-learn | 0.18 |
| scipy | 0.18.1 |
| seaborn | 0.7.0 |
| termcolor | 1.1.0 |

## Install

To install `scikit-ribo`, simply use the below command
    
    pip install scikit-ribo

## Usage

See the documentation on Read the Docs: [![Documentation Status](https://readthedocs.org/projects/scikit-ribo/badge/?version=latest)](http://scikit-ribo.readthedocs.io/en/latest/?badge=latest) or click [me](http://scikit-ribo.readthedocs.io/en/latest/)

For more information, please refer to the [template shell script](https://github.com/hanfang/scikit-ribo/blob/master/test/run_scikit_ribo.sh) about details of executing the two modules.

## Introduction

Scikit-ribo has two major modules: Ribosome A-site location prediction, and translation efficiency (TE) inference using a penalized generalized linear model (GLM). 

A complete analysis with scikit-ribo has two major procedures: 
1) data pre-processing to prepare the ORFs, codons for a genome: `scikit-ribo-build.py`
2) the actual model training and fitting: `scikit-ribo-run.py`

Inputs:
1) The alignment of Riboseq reads (bam)
2) Gene-level quantification of RNA-seq reads (from either Salmon or Kallisto)
3) A gene annotation file (gtf)
4) A reference genome for the model organism of interest (fasta)

Outpus:
1) Translation efficiency estimates for the genes
2) Translation elongation rate for 61 sense codons
3) Ribosome profile plots for each gene
4) Diagnostic plots of the models

## Reference

Fang et al, "Scikit-ribo: Accurate inference and robust modelling of translation dynamics at codon resolution" (Preprint coming up)
