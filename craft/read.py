import glob
import os

import pandas as pd
import numpy as np

def snptest(file):
    """ Read snptest data into an internal dataframe. """
    cols = ['chromosome','alleleA','alleleB','rsid','position','all_total', 'cases_total','controls_total','all_maf','frequentist_add_pvalue',
    'frequentist_add_beta_1', 'frequentist_add_se_1']
    df = pd.read_table(file, sep=' ', comment='#')[cols]
    df.rename(columns={'all_maf':'maf','frequentist_add_pvalue':'pvalue', 'frequentist_add_beta_1':'beta', 'frequentist_add_se_1':'se','alleleA':'allele1','alleleB':'allele2'}, inplace=True)
    return df

def plink(file, frq_file):
    """ Read plink (.assoc.logistic) data into an internal dataframe. """
    # read .assoc.logistic file
    cols = ['CHR','A1','SNP','BP','P','SE','OR']
    df = pd.read_csv(file, sep='\s+')[cols]
    df.rename(columns={'CHR':'chromosome','SNP':'rsid','BP':'position','A1':'allele1','P':'pvalue','SE':'se'}, inplace=True)
    # For finemap, we need the beta coefficient. For a binary logistic regression, ln(OR) = beta coefficient.
    for index, row in df.iterrows():
        df['beta'] = np.log(df['OR'])
    df.drop(columns='OR')

    # read .frq.cc file
    cols = ['CHR','SNP','A2','MAF_A','MAF_U','NCHROBS_A', 'NCHROBS_U']
    frq_df = pd.read_csv(frq_file, sep='\s+')[cols]
    frq_df.rename(columns={'CHR':'chromosome','SNP':'rsid','A2':'allele2','MAF_A':'maf','NCHROBS_A':'cases_total','NCHROBS_U':'controls_total'}, inplace=True)
    # if chromosome column has more than 1 number, read chromosome number from .assoc.logistic file and only include rows with that value.
    chromosomes = df.chromosome.unique()
    frq_df = frq_df[frq_df['chromosome'].isin(chromosomes)]
    frq_df = frq_df.drop("chromosome", axis=1)
    # create an all_total column
    for index, row in frq_df.iterrows():
        frq_df['all_total'] = frq_df['cases_total'] + frq_df['controls_total']
    # merge based on rsid
    df = pd.merge(df, frq_df, how='inner',on='rsid')
    # Rearrange column order after merge to match SNPtest format
    order = ['chromosome','allele1','allele2','rsid','position','all_total', 'cases_total','controls_total','maf','pvalue', 'beta', 'se']
    df = df[order]
    return df

def csv(file):
    """Read csv data into an internal dataframe. """
    cols = ['chromosome','rsid','alleleA','alleleB','position','all_total', 'cases_total','controls_total','all_maf','pvalue',
    'beta', 'se']
    df = pd.read_csv(file, sep='\t')[cols]
    return df

def maps(source_dir):
    """ Read genetic map data into a maps object. """
    map_file_list = glob.glob(source_dir + '/*chr[0-9]*.txt')
    maps = {}
    for file in map_file_list:
        map_file = pd.read_table(file)
        chromosome = map_file['Chromosome'].ix[0].strip('chr')
        maps[chromosome] = map_file
    return maps

def annovar(file , file_exonic, colnames):
    """ Read annovar output into an internal dataframe. """
    df = pd.DataFrame(columns=colnames)
    if os.path.getsize(file) != 0:
        df = df.append(pd.read_csv(file, sep='\t', names = colnames))
    if os.path.getsize(file_exonic) != 0:
        df2 = pd.read_csv(file_exonic, sep='\t', names = colnames, usecols=range(1,(len(colnames) + 1)))
        df = df.append(df2)
    return df
