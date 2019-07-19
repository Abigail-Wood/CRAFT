CRAFT (README)
==============

.. image:: http://readthedocs.org/projects/craft/badge/?version=latest
        :target: https://craft.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

Credible Refinement and Annotation of Functional Targets (CRAFT)

* Free, open-source software: MIT license
* Documentation: https://craft.readthedocs.io

This repository is still in development; our package is not due for full release until end July 2019.

(CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. It takes input as p-values from GWAS results. We have implemented this as a Python library, available via PyPI.

Quick start guide
-----------------
1. We recommend creating a virtual environment for package installation (craft and dependencies), using ``venv`` or conda.
2. Install Python 3.6 or later. Alternatively, you can create a local clone of this GitHub repository and run ``setup.py`` in your terminal.
3. Install python package at terminal using: ``python -m pip install --index-url https://test.pypi.org/simple/ --no-deps bio-craft``
4. Install ANNOVAR (and Perl if required, which ANNOVAR requires to run). http://annovar.openbioinformatics.org/en/latest/. You will also need to move the ANNOVAR directory into CRAFT, then download the hg19 database using the shell command ``annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGen -buildver hg19 humandb/``.
5. Optional - Download HapMap recombination maps for the correct human genome build from the NCBI FTP. Otherwise, you can use the version distributed by default with CRAFT (GRCh37/hg19, see genetic_maps directory for more information.)
6. Optional - Install supported finemapping packages you'd like to run (FINEMAP, PAINTOR or CAVIARBF) from their source websites and repositories (www.christianbrenner.com, https://github.com/gkichaev/PAINTOR_V3.0). If none are installed, you can use the approximate Bayes factor (ABF) calculation and credible SNP selection included in CRAFT.
7. Optional - Change the config file to match the locations and name of the ANNOVAR, finemapping packages and genetic maps directories (if required).
8. Optional - Make an empty 'output' folder if you want to run the test scripts.

The CRAFT workflow
------------------------------
.. image:: https://user-images.githubusercontent.com/15981287/60351501-d56bef00-99bd-11e9-8c6f-b4b6217c7d1b.png
        :alt: CRAFT workflow

Input file formats
------------------

GWAS summary statistics

| Input file is either in SNPTEST format or Plink.assoc.logistic format with an accompanying .frq.cc file.  
| Alternatively, you can reformat another dataset into the following columns and input it as a CSV file:  

| chromosome allele1 allele2 rsid position all_total cases_total controls_total maf pvalue beta se
| 2 T G rs76632663 43568820 12000 8415 3585 0.481 0.00000005 0.110147 0.18791
| 4 G A rs6823274 26098057 12000 8415 3585 0.345 0.00281 0.227728 0.155518 
| 6 T C rs28719598 128145388 11885 8300 3585 0.215 0.0000006 0.195838 0.139291
| 8 T C rs78519860 10943884 11885 8300 3585 0.101 0.5189 -0.154433 0.0791057


Output
------
Documentation still to be added; see output folder for examples of different types of output.

| CRAFT's ABF: produces an .abf.cred file as default.
| FINEMAP: produces .cred, .cred.annotated, .ld, .log_sss, .snp and .txt files as default.

Test data
---------
Currently we include two sets of test data distributed with the package, both for binary traits:

1. SNPTEST summary statistics from psoriatic arthritis (PsA) patients (chromosome 1)
2. PLINK .assoc.logistic and .frq.cc summary statistic files from PsA patients typed using Immunochip (chromosome 1)  

We have not (yet) tested this pipeline using data for quantitative traits, but have applied it to large datasets (>12 million SNPs) in patients with PsA, and are in the process of applying it in patients with JIA. 

Did you find an issue / missing feature?
----------------------------------------

We welcome all bug reports and requests for additional features using our GitHub issues tracker. This software will be supported by an active developer until at least 2020.

References
------------

Approximate Bayes Factor (ABF)

CRAFT uses a Python reimplemented version of the abf.R function written by [Chris Wallace](http://chr1swallace.github.io/) for the calculation of credible SNPs.

1. Jon Wakefield (2008) Bayes factors for genome-wide association studies: comparison with P-values. Genet Epidemiol DOI: 10.1002/gepi.20359
2. Bowes et. al (2015) Dense genotyping of immune-related susceptibility loci reveals new insights into the genetics of psoriatic arthritis.

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
