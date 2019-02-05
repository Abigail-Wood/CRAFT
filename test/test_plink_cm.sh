#!/bin/bash

python craft.py --file "test/plink_data/chr1.plink" --file_type plink --alpha 5e-5 --define_region cm --distance 0.1 --out plink_cm.test.craft