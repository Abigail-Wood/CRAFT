#!/bin/bash

python -m craft --file "test/psa_ichp_test/Immunochip_PsA_PhaseIII_gencall_QC_pca_corr_hg19_PC1-PC2_chr1.assoc.logistic" --type plink --frq "test/plink_data/basic.frq.cc" --alpha 5e-5 --distance_unit cm --distance 0.1 --outdir output
