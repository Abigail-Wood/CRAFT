#!/bin/bash

python -m craft --file "test/plink_data/basic.assoc.logistic" --type plink --frq "test/plink_data/basic.frq.cc" --alpha 5e-5 --distance_unit cm --distance 0.1 --out output/test.plink_cm_immunochip_ch1.index --outsf output/test.plink_cm_immunochip_ch1.cred
