Requirement
###########

Environment
-----------

- Python3
- Linux
- Recommend setting up your environment with `Conda <https://conda.io/docs/index.html>`_

Dependencies
------------

- Command-line pacakges:

+----------------+------------+
| Python package | Version >= |
+================+============+
| bedtools       | 2.26.0     |
+----------------+------------+

- Python package:

+----------------+------------+
| Python package | Version >= |
+================+============+
| colorama       | 0.3.7      |
+----------------+------------+
| glmnet_py      |0.1.0b      |
+----------------+------------+
| gffutils       | 0.8.7.1    |
+----------------+------------+
| matplotlib     | 1.5.1      |
+----------------+------------+
| numpy          | 1.11.2     |
+----------------+------------+
| pandas         | 0.19.2     |
+----------------+------------+
| pybedtools     | 0.7.8      |
+----------------+------------+
| pyfiglet       | 0.7.5      |
+----------------+------------+
| pysam          | 0.9.1.4    |
+----------------+------------+
| scikit_learn   | 0.18       |
+----------------+------------+
| scipy          | 0.18.1     |
+----------------+------------+
| seaborn        | 0.7.0      |
+----------------+------------+
| termcolor      | 1.1.0      |
+----------------+------------+

Note: When using pip install scikit-ribo, all the following dependencies will be pulled and installed automatically.

Installation
############

Options
-------
There are three options to install Scikit-ribo.


1. Install Scikit-ribo with pip::

    pip install scikit-ribo

2. Install Scikit-ribo with conda/biocodon::

    Coming up

3. Compile from source::

    git clone https://github.com/hanfang/scikit-ribo.git
    cd scikit-ribo
    python setup.py install

Test whether the installation is successful
-------------------------------------------
Once the installation is successful, you should expect the below if you type::

    scikit-ribo-run.py

.. image:: /images/successful_installation.png
   :align: center
   :scale: 75%