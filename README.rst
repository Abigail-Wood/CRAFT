CRAFT
================================================================

.. image:: http://readthedocs.org/projects/craft/badge/?version=latest
        :target: https://craft.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

Credible Refinement and Annotation of Functional Targets (CRAFT)

* Free, open-source software: MIT license
* Documentation: https://craft.readthedocs.io

PLEASE NOTE: 16/03/2019
=======================
This repository is still in development; our package is not due for release until May 2019.

(CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. It takes input as p-values from GWAS results. We have implemented this as a Python library, available via PyPI.

Quick start guide
-----------------
1. We recommend creating a virtual environment for package installation (craft and dependencies), using ``venv`` or conda.
2. Install Python 3.6 or later. Alternatively, you can create a local clone of the CRAFT GitHub repository and run ``setup.py`` in your terminal.
3. Install python package at terminal using: ``python -m pip install --index-url https://test.pypi.org/simple/ --no-deps bio-craft``
4. Install ANNOVAR_ (and Perl if required, which ANNOVAR requires to run)
.. _ANNOVAR: http://annovar.openbioinformatics.org/en/latest/
You will also need to move the ANNOVAR directory into CRAFT, then download the hg19 database using the shell command:
``annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGen -buildver hg19 humandb/ ``
5. Download HapMap recombination maps for the correct human genome build (e.g. GRCh37) from the NCBI FTP (or use the version distributed by default with CRAFT).
6. Install supported finemapping packages you'd like to run (FINEMAP or GCTA-COJO). If none are installed, you can use the default ABF calculation and credible SNP selection.
7. Optional: Change the config file to reflect the locations of the annovar, finemap and genetic maps directories if required.
8. Optional: Make an empty 'output' folder if you want to run the test scripts.

References
------------

Test data
---------
Currently two sets of test data: snptest format (chromosome 1, X SNPs) and PLINK format (chromosome 1, X SNPs) based on associations with psoriatic arthritis (PsA).

Input file formats
-----------------------

GWAS summary statistics

Input file is tab, space or comma-delimited and defined as follows:

``
SNPID      CHROM  POS       A1  A2  A1_UNAFF  PVAL
rs6823274   4     26098057  G   A   0.1484    0.4064
rs76632663  2     43568820  T   G   0.06988   0.4427
rs28719598  8     10943884  T   C   0.194     0.7702
rs78519860  6     128145388 T   C   0.07869   0.9007
``

Output
------

Did you find an issue / missing feature?
----------------------------------------

We welcome all bug reports and requests for additional features using our GitHub issues tracker. This software will be supported by an active developer until at least 2020.

Acknowledgements
----------------

CRAFT will use a reimplemented version of the abf.R function written by [Chris Wallace](http://chr1swallace.github.io/) for the calculation of credible SNPs.
