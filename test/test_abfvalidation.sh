#!/bin/bash

python craft.py --file "test/abf_validation_data/input.csv.csv" --file_type snptest --alpha 5e-5 --define_region cm --distance 0.1 --out snptest_cm.test.craft --outsf snptest_cm_abf.test.craft
