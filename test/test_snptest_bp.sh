#!/bin/bash

#make own output directory first in CRAFT folder
python -m craft --file "test/snptest_data/chr1.snptest.maf0.01.out" --type snptest --alpha 5e-5 --distance_unit bp --distance 500000 --out output/test.snptest_bp.index --outsf output/test.snptest_bp.cred
