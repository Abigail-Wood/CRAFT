#!/bin/bash

python -m craft --file "test/plink_data/chr1.plink" --type plink2 --alpha 5e-5 --distance_unit cm --distance 0.1 --out output/test.plink_cm.index --outsf output/test.plink_cm.cred
