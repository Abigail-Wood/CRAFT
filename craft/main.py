#!/usr/bin/env python
#
# Main CRAFT program

import sys
import os
import argparse
import glob

import pandas as pd

from craft import abf
from craft import annotate
from craft import config
from craft import log
from craft import read
from craft import finemap
from craft import paintor
from craft import visualise
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
        '--frq', action='store', help='Specify .frq file location (required for plink)')
    parser.add_argument(
        '--outdir', required=True, default='output',
        help='Specify output directory for all results - index set, credible set, ld file, finemapping output. Default = %(default)s.')
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
        help='For use with ABF, choose the cut-off threshold for cumulative posterior probability when determining credible sets. Default = %(default)s.')
    parser.add_argument(
        '--finemap_tool', choices={'finemap', 'paintor'},
        help='Choose which finemapping tool is used. Default = %(default)s.')
    parser.add_argument(
        '--n_causal_snps', type=int,
        help='For use with FINEMAP, specify the maximum number of causal snps considered in modelling. Default (set by FINEMAP) = 5')
    return parser.parse_args()

def main():
    options = parse_args() # Define command-line specified options
    file_names = glob.glob(options.file)
    if not file_names:
        log.error('Error: file not found!')
    for file in file_names:
        file_name = os.path.basename(os.path.normpath(file))
        file_dir = f"{options.outdir}/{file_name}"
        if os.path.exists(file_dir) == False:
            os.mkdir(file_dir)
        # Read input summary statistics
        if options.type == 'plink':
            if not options.frq:
                log.error('Error: .frq.cc file not found!')
            stats = read.plink(file, options.frq)
        else:
            reader = readers[options.type]
            stats = reader(file)
        # Get index SNPs
        if options.distance_unit == 'cm': # using cM as a distance unit
            distance = float(options.distance)
            maps = read.maps(config.genetic_map_dir)
            index_df = [gs.get_index_snps_cm(stats, options.alpha, distance, options.mhc, maps)]
        if options.distance_unit == 'bp': # using bp as a distance unit
            distance = int(options.distance)
            index_df = [gs.get_index_snps_bp(stats, options.alpha, distance, options.mhc)]
        index_df = pd.concat(index_df)

        # Output index SNPs. Float format is NOT default behaviour as this rounds to 6/7sf, use %g instead.
        index_df.to_csv(f"{os.path.join(file_dir, file_name)}.index", sep='\t', index=False, float_format='%g')

        # Get locus SNPs
        locus_dfs = gs.get_locus_snps(stats, index_df, options.distance_unit)

        # Calculate ABF and posterior probabilities
        data_list = abf.abf(locus_dfs, options.cred_threshold)

        # Annotate credible SNP set
        for data in data_list:
            data = annotate.prepare_df_annoVar(data)
            data = annotate.annotation_annoVar(data)
            # Output credible SNP set. Float format is NOT default behaviour as this rounds to 6/7sf, use %g instead.
            data.to_csv(f"{os.path.join(file_dir, data.index_rsid.unique()[0])}.abf.cred", sep='\t', index=False, float_format='%g')

        # Finemapping, if specified on command-line.
        if options.finemap_tool == "finemap":
            finemap.finemap(locus_dfs, index_df, file_dir, options.n_causal_snps)
            # Annotate finemap cred file results by iterating through index_df to find each .cred file
            i = 0
            for i, row in index_df.iterrows():
                cred_file = os.path.join(file_dir, row.rsid + ".cred")
                cred_dfs = read.finemap_cred(cred_file)
                cred_snps = pd.concat(cred_dfs)
                cred_snps_annotated = annotate.finemap_annotation_annoVar(cred_snps, locus_dfs[i])
                # write annotated SNPS dataframe as output file. Float format is NOT default behaviour as this rounds to 6/7sf, use %g instead.
                cred_snps_annotated.to_csv(f"{os.path.join(file_dir, row.rsid)}.cred.annotated", sep='\t', index=False, float_format='%g')
                # increment index count to select next index
                i+=1
        elif options.finemap_tool == "paintor":
            paintor.paintor(locus_dfs, index_df)

    return 0
