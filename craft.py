#!/usr/bin/env python

import sys
import os
import glob
import argparse
import pandas as pd
import numpy as np

# import craft module
import craft.craft_tools as cf
import craft.config as config

def get_options():

        parser = argparse.ArgumentParser()

        parser.add_argument('--file', action='store', dest='file', required=True,
                    help='Input summary statistics file. Can be multiple files with use of *')
        parser.add_argument('--bim', action='store', dest='bim_file',
                    help='PLINK bim file, required if using PLINK input')
        parser.add_argument('--snptest', action='store_const', const='snptest', dest='file_type',
                    help='Define input file type as snptest')
        parser.add_argument('--generic', action='store_const', const='generic', dest='file_type',
                    help='Define input file type as generic')
        parser.add_argument('--plink', action='store_const', const='plink', dest='file_type',
                    help='Define input file type as plink.assoc.logistic')
        parser.add_argument('--info', action='store', dest='info_file',
                    help='Tab delimited containing additional information to be merged with results (rsid required)')
        parser.add_argument('--alpha', action='store', dest='alpha', default=5e-8, type=float,
                    help='P-value threshold for declaring index SNPs. Default = 5e-8')
        parser.add_argument('--region_bp', action='store_const', const='bp', dest='region_type',
                    help='Define region around index SNP by base pairs')
        parser.add_argument('--region_cm', action='store_const', const='cm', dest='region_type',
                    help='Define region around index SNP by cM')
        parser.add_argument('--distance', action='store', dest='size', required=True,
                    help='Define distance around index SNP by base-pair or cM. E.g. 500000 for bp or 0.1 for cM')
        parser.add_argument('--no_mhc', action='store_false', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
        parser.add_argument('--mhc', action='store_true', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
        parser.add_argument('--out', action='store', dest='out_file', required=True,
                    help='Output path')
        parser.set_defaults(mhc=False)
        
        args =  parser.parse_args()

        # check bim has been specified with plink
        if args.file_type == 'plink' and not args.bim_file:
            parser.error('Input Error. A bim file is required when using plink input.')
        return args

def get_options_2():

        parser = argparse.ArgumentParser()

        parser.add_argument('--file', action='store', dest='file', required=True,
                    help='Input summary statistics file. Can be multiple files with use of *')
        parser.add_argument('--bim', action='store', dest='bim_file',
                    help='PLINK bim file, required if using PLINK input')
        parser.add_argument('--file_type', action='store', choices=['snptest','plink','generic'], dest='file_type',
                    help='Define input file type - snptest, plink, generic')
        parser.add_argument('--decorate', action='store', dest='decorate_file',
                    help='Tab delimited containing additional information to be merged with results (rsid required)')
        parser.add_argument('--alpha', action='store', dest='alpha', default=5e-8, type=float,
                    help='P-value threshold for declaring index SNPs. Default = 5e-8')
        parser.add_argument('--define_region', action='store', choices=['cm','bp'], dest='region_type',
                    help='Define region around index SNP cm or bp')
        parser.add_argument('--distance', action='store', dest='size', required=True,
                    help='Define distance around index SNP by base-pair or cM. E.g. 500000 for bp or 0.1 for cM')
        parser.add_argument('--no_mhc', action='store_false', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
        parser.add_argument('--mhc', action='store_true', dest='mhc',
                    help='Exclude the MHC region. Default is to exclude')
        parser.add_argument('--out', action='store', dest='out_file', required=True,
                    help='Output path')
        parser.set_defaults(mhc=False)
        
        args =  parser.parse_args()

        # check bim has been specified with plink
        #if args.file_type == 'plink' and not args.bim_file:
        #    parser.error('Input Error. A bim file is required when using plink input.')
        return args


def main():

        # define input/output
        options = get_options_2()

        # read summary stats
        file_names = glob.glob(options.file)

        if options.file_type == "snptest":
                stats_list = map(cf.read_snptest, file_names)

        if options.file_type == "generic":
                stats_list = map(cf.read_generic, file_names)

        if options.file_type == "plink":
                stats_list = map(cf.read_plink_2, file_names)

        # get index snps
        if options.region_type == "bp":
                 size = int(options.size)
                 index_list = map(lambda d : cf.get_index_snps_bp(d, options.alpha, size, options.mhc), stats_list)

        if options.region_type == "cm":
                 maps = cf.read_maps(config.genetic_map_dir)
                 size = float(options.size)
                 index_list = map(lambda d : cf.get_index_snps_cm(d, options.alpha, size, options.mhc, maps), stats_list)

        # create df of index snps
        index_df = pd.concat(index_list)

        print index_df.head()

        index_df.chromosome = index_df.chromosome.astype(int)
        index_df.position = index_df.position.astype(int)

        # annotate index SNPs
        annotated_index_df = cf.base_annotation(index_df)

        # add additional information if supplied
        if options.decorate_file:
            annotated_index_df = cf.merge_info(annotated_index_df, options.decorate_file)

        # write index SNP table
        out_file = options.out_file
        annotated_index_df.to_csv(out_file, sep='\t', index=False)

if __name__=='__main__':
	main()
