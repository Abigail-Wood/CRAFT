#!/usr/bin/env python
#
# Main CRAFT program

import sys
import os
import argparse
import glob

import pandas as pd

from craft import abf_beta
from craft import annotate
from craft import config
from craft import log
from craft import read
from craft import finemap
import craft.getSNPs as gs

# All file reading functions. Each takes a file name and returns a DataFrame.
readers = {'snptest': read.snptest,
           'plink': read.plink,
           'csv': read.csv,
}

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--file', required=True,
        help='Input summary statistics file. Use * to include multiple files (must have the same file type.)')
    parser.add_argument(
        '--type', required = True, choices=readers.keys(),
        help='Define input file type.')
    parser.add_argument(
        '--bim', action='store', help='Specify .bim file location (required for plink). Use * to include multiple files (must be in the same order as matching plink files)')
    parser.add_argument(
        '--frq', action='store', help='Specify .frq file location (required for plink)')
    parser.add_argument(
        '--out', required=True,
        help='Output file for table of index SNPs with summary statistics.')
    parser.add_argument(
        '--outsf', required=True,
        help='Output file for summary statistics of the credible set.')
    parser.add_argument(
        '--alpha', default=5e-8, type=float,
        help='P-value threshold for declaring index SNPs. Default = %(default)s.')
    parser.add_argument(
        '--distance_unit', choices=['cm','bp'], default = 'cm',
        help='Choose the distance unit (for use in defining regions). Default = %(default)s.')
    parser.add_argument(
        '--distance', default='0.1',
        help='Define distance around index SNP by base-pair or cM. e.g. 500000 for bp or 0.1 for cM. Default = %(default)s.')
    parser.add_argument(
        '--mhc', action='store_true',
        help='Include the MHC region. Default = %(default)s.')
    parser.add_argument(
        '--cred_threshold', choices={'95', '99'}, default='95',
        help='Choose the cut-off threshold for cumulative posterior probability, when determining credible sets with ABF. Default = %(default)s.')
    parser.add_argument(
        '--finemap_tool', choices={'cojo', 'finemap', 'paintor'},
        help='Choose which finemap tool is used. Default = %(default)s.')
    return parser.parse_args()

def main():
    options = parse_args() # Define command-line specified options
    file_names = glob.glob(options.file)
    if len(file_names) != 1:
        log.error("Can't (yet) process more than one file.")
    #TODO: extend to take multiple input files
    if not file_names:
        log.error('Error: file not found!')

    # Read input summary statistics
    if options.type == 'plink': # dirty hack
        stats = [read.plink(n, options.frq) for n in file_names]
    else:
        reader = readers[options.type]
        stats = [reader(n) for n in file_names]

    # Get index SNPs
    if options.distance_unit == 'cm': # using cM as a distance unit
        distance = float(options.distance)
        maps = read.maps(config.genetic_map_dir)
        index_dfs = [gs.get_index_snps_cm(d, options.alpha, distance, options.mhc, maps) for d in stats]
    if options.distance_unit == 'bp': # using bp as a distance unit
        distance = int(options.distance)
        index_dfs = [gs.get_index_snps_bp(d, options.alpha, distance, options.mhc) for d in stats]
    index_df = pd.concat(index_dfs)
    out_file = options.out

    # Output index SNPs
    index_df.to_csv(out_file, sep='\t', float_format='%5f', index=False)

    # Get locus SNPs
    for stat_df in stats:
        data_dfs = gs.get_locus_snps(stat_df, index_df, options.distance_unit)

    # Calculate ABF and posterior probabilities
    data_list = abf_beta.abf(data_dfs, options.cred_threshold)
    data = pd.concat(data_list)

    # Annotate credible SNP set
    data = annotate.prepare_df_annoVar(data)
    data = annotate.base_annotation_annoVar(data) # Annotate credible SNPs

    # Output credible SNP set
    data.to_csv(options.outsf, sep='\t', float_format='%5f', index=False)

    # Finemapping, if specified on command-line.
    if options.finemap_tool == "finemap":
        finemap.finemap(data_dfs, index_df)
    elif options.finemap_tool == "paintor":
        paintor.paintor(data_dfs, index_df)
    return 0
