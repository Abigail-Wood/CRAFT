#!/bin/bash

python -m craft --file "test/PLINK_test/Immunochip_PsA_ch1_PLINK.assoc.logistic" --type plink --frq frq_file.frq --alpha 5e-5 --distance_unit cm --distance 0.1 --out output/test.plink_cm_immunochip_ch1.index --outsf output/test.plink_cm_immunochip_ch1.cred
