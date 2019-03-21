#!/bin/bash

python -m craft --file "test/snptest_data/chr1.snptest.maf0.01.out" --type snptest --alpha 5e-5 --distance_unit bp --distance 500000 --out test.snptest_bp.index --outsf test.snptest_bp.cred
