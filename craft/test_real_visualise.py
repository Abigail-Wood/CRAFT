#!/usr/bin/env python

import sys
import argparse
import pdb

import scipy.stats
from pandas import DataFrame

import numpy

import visualise
import read

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', '-d', action='store_true',
                        help='If any error arises, fall into the Python debugger.')
    parser.add_argument('--verbose', '-v', action='count',
                        help='Increase logging verbosity..')
    parser.add_argument('--input_file', '-i',
                        help='Locate original input file.')
    parser.add_argument('--index_file', '-ix',
                        help='Locate CRAFT .index output file')
    parser.add_argument('--test', '-t', action='store_true',
                        help='Draw test charts.')
    parser.add_argument('--ld', '-l', action='store_true',
                        help='Draw LD block charts.')
    parser.add_argument('--manhattan', '-m', action='store_true',
                        help='Draw Manhattan charts.')
    parser.add_argument('--locus', '-c', action='store_true',
                        help='Draw locus charts.')

    return parser.parse_args()

def run(options):

    if options.manhattan:
        # read in dataframe
        df = read.snptest(options.input_file)
        index_df = read.index(options.index_file)
        # using all the default parameters.
        manhattan_chart1 = visualise.manhattan(df, f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5, index_df=index_df)
        manhattan_chart1.savefig('manhattan1.png', dpi=300)
        # no alpha: no SNPs are marked as special.
        manhattan_chart2 = visualise.manhattan(df, f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5, index_df=index_df)
        manhattan_chart2.savefig('manhattan2.png', dpi=300)
        # draw a load of things differently. Good SNPs are distinguished
        # but not labelled.
        manhattan_chart3 = visualise.manhattan(df, f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5, index_df=index_df,
        color='g', good_color='r',
        size=0.5, good_size=25,
        marker='o', good_marker='x',
        good_label_column=None)
        manhattan_chart3.savefig('manhattan3.png', dpi=300)
        # rotating labels and skipping the alpha line
        manhattan_chart4 = visualise.manhattan(df, f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5, index_df=index_df,
        alpha_line_style='',
        good_label_rotation=45)
        manhattan_chart4.savefig('manhattan4.png', dpi=300)

def main():
    options = parse_args() # Define command-line specified options
    try:
        run(options)
    except:
        if options.debug:
            pdb.post_mortem()
        else:
            raise

if __name__=='__main__':
	main()
