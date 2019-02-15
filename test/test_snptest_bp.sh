#!/bin/bash

python craft.py --file "test/snptest_data/chr1.snptest.maf0.01.out" --file_type snptest --alpha 5e-6 --define_region bp --distance 500000 --out snptest_bp.test.craft
