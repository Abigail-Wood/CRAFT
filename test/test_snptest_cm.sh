#!/bin/bash

python craft.py --file "test/snptest_data/chr1.snptest.maf0.01.out" --file_type snptest --alpha 5e-5 --define_region cm --distance 0.1 --out snptest_cm.test.craft