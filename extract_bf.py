#!/usr/bin/env python3

# Goal: extract all lines that read '- Log10-BF of >= one causal SNP : 0.476 in [rsid].log_sss files, then match rsid in filename to index line in .index file and append the value in a new column at the end of the row.

import argparse
import sys
import glob

import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input',required=True,nargs='*',type=str,
                        help='List of log_sss files.')
    parser.add_argument('-ix','--index',required=True,type=str,
                        help='Index file name.')
    parser.add_argument('-o','--output',nargs='?',default=sys.stdout, type=str,
                        help='Output name for index file')
    return parser.parse_args()

def extract_bf(filenames):
    start = 'out/'
    end = '.'
    index_bf = []
    for file in filenames:
        index = file[file.find(start)+len(start):file.rfind(end)]
        with open(file,'r') as f:
            for line in f:
                if line.startswith('- Log10-BF'):
                    bf = line.split(':')[1]
                    bf = bf.strip()
                    index_bf.append((index,bf))
    return index_bf

def add_bf_to_index(file, index_bf,output):
    index_df = pd.read_csv(file, sep='\s+')
    index_df.insert(len(index_df.columns),'Log10-BF', 0)
    for index, bf in index_bf:
        index_df.loc[(index_df['rsid'] == index), 'Log10-BF'] = bf
    index_df.to_csv(f'{output}.csv', sep='\t',index=False)

def main():
    args = parse_args()
    index_bf = extract_bf(args.input)
    add_bf_to_index(args.index, index_bf, args.output)

if __name__ == "__main__":
    main()
