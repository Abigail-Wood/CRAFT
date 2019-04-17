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
                        help='''If any error arises, fall into the Python
                        debugger.''')
    parser.add_argument('--verbose', '-v', action='count',
                        help='Increase logging verbosity..')
    parser.add_argument('--input_file', '-i',
                        help='Locate original input file.')
    parser.add_argument('--index_file', '-ix',
                        help='Locate CRAFT .index output file.')
    parser.add_argument('--ld_file', '-ldi',
                        help='Locate LD input file for locus.')
    parser.add_argument('--ld_rsids', '-ldr',
                        help='Locate list of RSIDs (variant_file.txt)')
    parser.add_argument('--cred_file', '-c',
                        help='Locate .cred file (abf or finemap.)')
    parser.add_argument('--manhattan', '-m', action='store_true',
                        help='Draw Manhattan charts.')
    parser.add_argument('--ld', '-l', action='store_true',
                        help='Draw LD block charts.')
    parser.add_argument('--locus', '-u', action='store_true',
                        help='Draw locus charts.')
    parser.add_argument('--cred_type', '-t', choices=['finemap','abf'],
                        help='Specify type of .cred file provided.')

    return parser.parse_args()

def run(options):
    if options.manhattan:
        # read in dataframe
        df = read.snptest(options.input_file)
        index_df = read.index(options.index_file)
        # display in green with red x, no rsid labels.
        manhattan_green = visualise.manhattan(df,
        f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5,
        index_df=index_df, color='g', good_color='r', size=0.5, good_size=25,
        marker='o', good_marker='x', good_label_column=None)
        manhattan_chart3.savefig("manhattan__1.png", dpi=300)

        # rotating labels and skipping the alpha line
        manhattan_blue = visualise.manhattan(df,
        f"chromosome {df.chromosome.unique()[0]}", alpha=5e-5,
        index_df=index_df, alpha_line_style='', good_label_rotation=45)
        manhattan_chart4.savefig("manhattan_2.png", dpi=300)

    if options.ld:
        ld_array = read.ld(options.ld_file)
        cred_df = read.cred_annotated(options.cred_file)
        cred_df = cred_df.sort_values('position', ascending=True)
        cred_snps = list(cred_df['rsid'])
        # build dictionary of SNP name to index from variant file
        variant_df = read.variant_file(options.ld_rsids)
        variant_dict = dict(zip(variant_df['RSID'],variant_df.index))
        indexes = []
        for key in variant_dict:
            indexes.append(variant_dict[key])
        positions = variant_df.position.unique()

        ld_chart1 = visualise.ld_block(ld_array, indexes, names=None,
                    labels=dict(mid=f"chromosome {cred_df.chromosome.unique()[0]}", left=min(positions), right=max(positions)))
        ld_chart1.savefig('ld1.png', dpi=300)

        indexes = []
        for snp in cred_snps:
            assert snp in variant_dict
            indexes.append(variant_dict[snp])
        # add positions and names using cred SNPs annotated.
        positions = cred_df.position.unique()
        names = cred_snps

        ld_chart2 = visualise.ld_block(ld_array, indexes, names,
                    labels=dict(mid=f"chromosome {cred_df.chromosome.unique()[0]}", left=min(positions), right=max(positions)))
        ld_chart2.savefig('ld2.png', dpi=300)

    if options.locus and options.cred_type == 'abf':
        # tracknames come from unique values in var_effect column
        cred_df = read.abf_cred(options.cred_file)
        tracknames= cred_df.var_effect.unique()
        # Don't yet understand this bit.
        nt = scipy.stats.hypergeom(50, len(tracknames), len(tracknames))
        tracks = [','.join(numpy.random.choice(tracknames,size=i,replace=False)) for i in nt.rvs(size=1196)]
        # Dataframe for cred snps
        df = DataFrame({'rsid' : cred_df['rsid'],
                        'position' : cred_df['position'],
                        'posterior' : cred_df['pp'],
                        'tracks' : tracks})
        # with the probability chart on the bottom
        locus_chart1 = visualise.locus(df, tracks=None, track_column='tracks', pos_top = False)
        locus_chart1.savefig('locus1.png', dpi=300)
        # without any tracks chart
        locus_chart2 = visualise.locus(df, tracks=False)
        locus_chart2.savefig('locus2.png', dpi=300)
        # With the probability chart on the bottom and a number of other parameter tweaks.
        # This shows some of the other design possibilities.
        locus_chart3 = visualise.locus(df,
                                    alpha_line_color='g',
                                    alpha_line_style=':',
                                    good_size=10,
                                    good_marker='x',
                                    good_color='r',
                                    good_label_column=None,
                                    tracks=tracknames[2:5],
                                    track_good_linelength=0.8,
                                    track_linelength=0.8,
                                    track_alpha = 1,
                                    pos_top = False,
                                    track_height=0.25,
                                    figsize=(6,6),
                                    track_lines = True)
        locus_chart3.savefig('locus3.png', dpi=300)

    if options.locus and options.cred_type =='finemap':
        # tracknames come from unique values in var_effect column
        cred_df = read.finemap_cred(options.cred_file)
        tracknames= cred_df.var_effect.unique()
        # Don't yet understand this bit.
        nt = scipy.stats.hypergeom(50, len(tracknames), len(tracknames))
        tracks = [','.join(numpy.random.choice(tracknames,size=i,replace=False)) for i in nt.rvs(size=1196)]
        # Dataframe for cred snps
        df = DataFrame({'rsid' : cred_df['rsid'],
                        'position' : cred_df['position'],
                        'posterior' : cred_df['pp'],
                        'tracks' : tracks})
        # with the probability chart on the bottom
        locus_chart1 = visualise.locus(df, tracks=None, track_column='tracks', pos_top = False)
        locus_chart1.savefig('locus1.png', dpi=300)
        # without any tracks chart
        locus_chart2 = visualise.locus(df, tracks=False)
        locus_chart2.savefig('locus2.png', dpi=300)
        # With the probability chart on the bottom and a number of other parameter tweaks.
        # This shows some of the other design possibilities.
        locus_chart3 = visualise.locus(df,
                                    alpha_line_color='g',
                                    alpha_line_style=':',
                                    good_size=10,
                                    good_marker='x',
                                    good_color='r',
                                    good_label_column=None,
                                    tracks=tracknames[2:5],
                                    track_good_linelength=0.8,
                                    track_linelength=0.8,
                                    track_alpha = 1,
                                    pos_top = False,
                                    track_height=0.25,
                                    figsize=(6,6),
                                    track_lines = True)
        locus_chart3.savefig('locus3.png', dpi=300)

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
