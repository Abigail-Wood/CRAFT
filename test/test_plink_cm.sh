#!/bin/bash

python craft.py --file "test/plink_data/chr1.plink" --type plink2 --alpha 5e-5 --distance_unit cm --distance 0.1 --out test.plink_cm.index --outsf test.plink_cm.cred
