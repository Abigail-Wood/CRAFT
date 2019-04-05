#!/bin/bash

python -m craft --file "test/snptest_data/chr1.snptest.maf0.01.out" --type snptest --alpha 5e-5 --distance_unit cm --distance 0.1 --out output/test.snptest_cm_paintor.index --outsf output/test.snptest_cm_paintor.cred --finemap_tool paintor
