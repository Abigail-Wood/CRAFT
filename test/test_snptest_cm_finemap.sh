#!/bin/bash

python -m craft --file "test/snptest_data/chr1.snptest.maf0.01.out" --type snptest --alpha 5e-5 --distance_unit cm --distance 0.1 --outdir output/ --finemap_tool finemap --n_causal_snps 3
