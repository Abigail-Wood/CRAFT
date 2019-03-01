#!/usr/bin/env python

import sys
import os
import argparse
import glob

import pandas as pd
import numpy as np

import craft.config as config
import craft.read as read
import craft.getSNPs as gs
import craft.abf as abf
import craft.annotate as annotate
import craft.finemap as finemap

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
    parser.add_argument('--finemap', action='store_true', dest='finemap',
                    help='Output file path for sum stats with credible set')
    parser.add_argument('--no-finemap', action='store_false', dest='finemap',
                    help='Output file path for sum stats with credible set')
    parser.set_defaults(mhc=False, finemap=True)

    # TO-DO: Add argument testing, e.g. check bim specified with plink
    return parser.parse_args()

# All the reading functions. Each takes a file name and returns a DataFrame.

readers = {'snptest': read.snptest,
           'plink': read.plink,
           'plink2': read.plink_noBIM,
           'indexsnps': read.indexsnps,
           'generic': read.generic,
}

def main():
    # Define all command-line specified options
    options = parse_args()

    # Read input summary statistics file (according to format)
    file_names = glob.glob(options.file)
    reader = readers[options.file_type]
    stats = [reader(n) for n in file_names]

    # Get index SNPs
    if options.region_type == 'cm':
        size = float(options.size)
        maps = read.maps(config.genetic_map_dir)
        index = [gs.get_index_snps_cm(d, options.alpha, size, options.mhc,     maps) for d in stats]
    if options.region_type == 'bp':
        size = int(options.size)
        index = [gs.get_index_snps_bp(d, options.alpha, size, options.mhc)
                for d in stats]

    index_df = pd.concat(index) # Create df of index SNPs
    index_df.chromosome = index_df.chromosome.astype(int)
    index_df.position = index_df.position.astype(int)

    # Annotate index SNPs (commented out)
    ## annotated_index_df = gs.base_annotation(index_df)
    # Add additional information (if supplied)
    # if options.decorate_file:
    #    annotated_index_df = annotate.merge_info(annotated_index_df, options.decorate_file)

    # Write index SNP table to specified output file
    out_file = options.out_file
    index_df.to_csv(out_file, sep='\t', index=False)

    # If stats analysis output filename, calculate ABF and posterior prob.
    if options.outsf:
        for stat_df in stats:
            data = gs.get_locus_snps(stat_df, index_df)
        data['ABF'] = data.apply(
            lambda row: abf.calc_abf(pval=row['pvalue'],
                                maf=row['all_maf'],
                                n=row['all_total'],
                                n_controls=row['controls_total'],
                                n_cases=row['cases_total']), axis=1)
        data = data.sort_values('ABF', ascending=False)
        data = abf.calc_postprob(data)
        data = abf.calc_postprobsum(data)

        # Write all data to output file
        data.to_csv(options.outsf, sep='\t', float_format='%5f', index=False)

    if options.finemap:
        print('finemap is on!')

    return 0

if __name__=='__main__':
    main()
