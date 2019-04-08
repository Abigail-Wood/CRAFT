CRAFT
================================================================

.. image:: http://readthedocs.org/projects/craft/badge/?version=latest
        :target: https://craft.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

Credible Refinement and Annotation of Functional Targets (CRAFT)

* Free, open-source software: MIT license
* Documentation: https://craft.readthedocs.io

PLEASE NOTE: 04/04/2019
=======================
This repository is still in development; our package is not due for release until May 2019.

(CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. It takes input as p-values from GWAS results. We have implemented this as a Python library, available via PyPI.

Quick start guide
-----------------
1. We recommend creating a virtual environment for package installation (craft and dependencies), using ``venv`` or conda.
2. Install Python 3.6 or later. Alternatively, you can create a local clone of this GitHub repository and run ``setup.py`` in your terminal.
3. Install python package at terminal using: ``python -m pip install --index-url https://test.pypi.org/simple/ --no-deps bio-craft``
4. Install ANNOVAR
5. (and Perl if required, which ANNOVAR requires to run). http://annovar.openbioinformatics.org/en/latest/. You will also need to move the ANNOVAR directory into CRAFT, then download the hg19 database using the shell command:
``annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGen -buildver hg19 humandb/``
5. Optional: Download HapMap recombination maps for the correct human genome build (e.g. GRCh37) from the NCBI FTP. Otherwise, you will use the version distributed with CRAFT (see genetic_maps directory for more information.)
6. Optional: Install supported finemapping packages you'd like to run (FINEMAP, PAINTOR or CAVIARBF) from their source websites and repositories (www.christianbrenner.com, https://github.com/gkichaev/PAINTOR_V3.0). If none are installed, you can still use by default the approximate Bayes factor (ABF) calculation and credible SNP selection.
7. Optional: Change the config file to match the locations and name of the ANNOVAR, finemapping packages and genetic maps directories (if required).
8. Optional: Make an empty 'output' folder if you want to run the test scripts.

Input file formats
------------------

GWAS summary statistics

Input file is tab, space or comma-delimited and defined as follows:

``
RSID      CHROM  POS       A1  A2  A1_UNAFF  PVAL
rs6823274   4     26098057  G   A   0.1484    0.4064
rs76632663  2     43568820  T   G   0.06988   0.4427
rs28719598  8     10943884  T   C   0.194     0.7702
rs78519860  6     128145388 T   C   0.07869   0.9007
``

Output
------

Test data
---------
Currently we include two sets of test data, both for binary traits:
1. SNPTEST summary statistics from PsA patients (chromosome 1)
2. PLINK .assoc.logistic and .frq.cc summary statistic files from PsA patients typed using Immunochip (chromosome 1)
We have not (yet) tested this pipeline using data for quantitative traits.

Did you find an issue / missing feature?
----------------------------------------

We welcome all bug reports and requests for additional features using our GitHub issues tracker. This software will be supported by an active developer until at least 2020.

References
------------

Approximate Bayes Factor (ABF)
CRAFT uses a Python reimplemented version of the abf.R function written by [Chris Wallace](http://chr1swallace.github.io/) for the calculation of credible SNPs.
1. Jon Wakefield (2008) Bayes factors for genome-wide association studies: comparison with P-values. Genet Epidemiol DOI: 10.1002/gepi.20359
2. http://www.nature.com/ncomms/2015/150205/ncomms7046/abs/ncomms7046.html](Bowes et. al (2015) Dense genotyping of immune-related susceptibility loci reveals new insights into the genetics of psoriatic arthritis).

ANNOVAR
1. Wang, K., Li, M. & Hakonarson, H. ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data. Nucleic Acids Res. 38, e164 (2010).

CAVIARBF (C++ version, v0.1.4.1)
1. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4512539/

FINEMAP(v1.3)

1. Refining fine-mapping: effect sizes and regional heritability. bioRxiv. (2018).
2. Prospects of fine-mapping trait-associated genomic regions by using summary statistics from genome-wide association studies. Am. J. Hum. Genet. (2017).
3. FINEMAP: Efficient variable selection using summary data from genome-wide association studies. Bioinformatics 32, 1493-1501 (2016).

PAINTOR (v3.0)
1. Kichaev et al. (PLOS Genetics, 2014)
2. Kichaev et al. (American Journal of Human Genetics, 2015)
3. Kichaev et al. (Bioinformatics, 2016).
