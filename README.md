CRAFT
=====

**C**redible **R**efinement and **A**nnotation of **F**unctional **T**argets (CRAFT) is a pipeline for the calculation, annotation and visualisation of credible SNP sets. It takes input as p-values from GWAS results.

We have implemented this as a Python library, with the following sets of modules:
* Input_data: modules for importing PLINK, SNPTEST data

Quick-start Guide
------------

1. Install python package at terminal using: `python install bio-craft`

Requirements
------------

__Software__ (TO BE VERIFIED)
* Python 3
    - vcf  
    - pandas  
    - numpy
    - scipy

Source data
-----------

__Variant Effect Predictor (VEP)__

Get the GRCh37 (release 84) core cache (not included due to size)
```
mkdir -p source_data/ensembl/{cache,plugins}
cd source_data/ensembl/cache/
wget ftp://ftp.ensembl.org/pub/release-84/variation/VEP/homo_sapiens_vep_84_GRCh37.tar.gz
tar -zxvf homo_sapiens_vep_84_GRCh37.tar.gz
cd -
```

__HapMap recombination map__

```
mkdir -p source_data/genetic_map_HapMapII_GRCh37/
cd source_data/genetic_map_HapMapII_GRCh37/
wget ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
tar -zxvf genetic_map_HapMapII_GRCh37.tar.gz
cd -
```

Test data
---------
Currently two sets of test data: snptest format (chromosome 1, X SNPs) and plink format (chromosome 1, X SNPs) based on associations with psoriatic arthritis (PsA).

User input file formats
-----------------------

__GWAS summary statistics__

Input file is tab, space or comma-delimited and defined as folllows:

```
SNPID      CHROM  POS       A1  A2  A1_UNAFF  PVAL
rs6823274   4     26098057  G   A   0.1484    0.4064
rs76632663  2     43568820  T   G   0.06988   0.4427
rs28719598  8     10943884  T   C   0.194     0.7702
rs78519860  6     128145388 T   C   0.07869   0.9007
```

Output
------

Did you find an issue / missing feature?
----

We welcome all bug reports and requests for additional features using our GitHub issues tracker. This software will be supported by an active developer until at least 2020.

Acknowledgements
----------------

CRAFT will use a reimplemented version of the abf.R function written by [Chris Wallace](http://chr1swallace.github.io/) for the calculation of credible SNPs.
