#!/usr/bin/env python

import sys
import os
import argparse
import glob

import pandas as pd
import numpy as np

import craft.config as config
import craft.read as read
import craft.getIndexSNPs as cf
import craft.abf as abf
import craft.annotate as annotate

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', action='store', dest='file', required=True,
                    help='Input summary statistics file. Can be multiple files with use of *')
    parser.add_argument('--file_type', action='store', choices=['snptest','plink','plink2','generic'], dest='file_type',
                    help='Define input file type - snptest, plink, plink2, indexsnps, generic')
    parser.add_argument('--bim', action='store', dest='bim_file',
                    help='Specify BIM file location (required if using plink input)')
    parser.add_argument('--decorate', action='store', dest='decorate_file',
                    help='Tab delimited file containing additional information to be merged with results (rsid required)')
    parser.add_argument('--alpha', action='store', dest='alpha', default=5e-8, type=float,
                    help='P-value threshold for declaring index SNPs. Default = 5e-8')
    parser.add_argument('--define_region', action='store', choices=['cm','bp'], dest='region_type',
                    help='Define region around index SNP cm or bp')
    parser.add_argument('--distance', action='store', dest='size', required=True,
                    help='Define distance around index SNP by base-pair or cM. e.g. 500000 for bp or 0.1 for cM')
    parser.add_argument('--no_mhc', action='store_false', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
    parser.add_argument('--mhc', action='store_true', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
    parser.add_argument('--out', action='store', dest='out_file', required=True,
                    help='Output file path for index SNP table')
    parser.add_argument('--outsf', metavar='<file>', type=str, dest='outsf',
                    help='Output file path for sum stats with credible set')
    parser.set_defaults(mhc=False)

    # TO-DO: Add argument testing, e.g. check bim specified with plink
    return parser.parse_args()

def main():
    # Define all command-line specified options
    options = parse_args()

    # Read input summary statistics file (according to format)
    file_names = glob.glob(options.file)
    if options.file_type == "snptest":
        stats_list = map(read.snptest, file_names)
    if options.file_type == "plink":
        stats_list = map(read.plink, file_names)
    if options.file_type == "plink_noBIM":
        stats_list = map(read.plink2, file_names)
    if options.file_type == "indexsnps":
        stats_list = map(read.indexsnps, file_names)
    if options.file_type == "generic":
        stats_list = map(read.generic, file_names)

    # Get index SNPs
    if options.region_type == "cm":
        maps = read.maps(config.genetic_map_dir)
        size = float(options.size)
        index_list = map(lambda d : cf.get_index_snps_cm(d, options.alpha, size, options.mhc, maps), stats_list)

    if options.region_type == "bp":
        size = int(options.size)
        index_list = map(lambda d : cf.get_index_snps_bp(d, options.alpha, size, options.mhc), stats_list)

    # Create df of index SNPs
    index_df = pd.concat(index_list)

    print(index_df.head())

    index_df.chromosome = index_df.chromosome.astype(int)
    index_df.position = index_df.position.astype(int)

    # Annotate index SNPs (commented out)
    ## annotated_index_df = cf.base_annotation(index_df)

    # Add additional information (if supplied)
    if options.decorate_file:
        annotated_index_df = annotate.merge_info(annotated_index_df, options.decorate_file)

    # Write index SNP table to specified output file
    out_file = options.out_file
    index_df.to_csv(out_file, sep='\t', index=False)

    # If stats analysis output file name specified then run ABF
    if options.outsf != None:
        data = index_df
        data["logABF"] = data.apply(
            lambda row: abf.calc_abf(pval=row['pvalue'],
                                 maf=row['all_maf'],
                                 n=row['all_total'],
                                 prop_cases=None), axis=1)
        data = data.sort_values("logABF", ascending=False)

        # Calculate posterior probability for each SNP
        sum_lABF = abf.log_sum(data["logABF"])
        data["postprob"] = data["logABF"].apply(np.exp) / np.exp(sum_lABF)

        # Calc cumulative sum of the posterior probabilities
        data["postprob_cumsum"] = data["postprob"].cumsum()

        # Find 99% and 95% credible sets - this is horrible
        set_idx = data["postprob_cumsum"].gt(0.99).tolist().index(True)
        data["is99_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))
        set_idx = data["postprob_cumsum"].gt(0.95).tolist().index(True)
        data["is95_credset"] = [1] * (set_idx + 1) + [0] * (data.shape[0] - (set_idx + 1))

        # Write
        data.to_csv(options.outsf, sep='\t', index=False)

        return 0

if __name__=='__main__':
    main()
