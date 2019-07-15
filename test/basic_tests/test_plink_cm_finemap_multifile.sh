#!/bin/bash

python -m craft --file "test/psa_ichp_test/*.assoc.logistic" --type plink --frq "test/psa_ichp_test/Immunochip_PsA_PhaseIII_gencall_QC_hg19_updateNames.frq.cc" --alpha 5e-5 --distance_unit cm --distance 0.1 --outdir output --finemap_tool finemap -n_causal_snps 3
