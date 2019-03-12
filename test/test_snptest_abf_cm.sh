#!/bin/bash

python craft.py --file "test/snptest_data/chr1.snptest.maf0.01.out" --file_type snptest --alpha 5e-5 --define_region cm --distance 0.1 --out test.snptest_cm.index --outsf test.snptest_cm_abf --no-finemap
