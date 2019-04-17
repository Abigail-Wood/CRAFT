import glob
import os

import pandas as pd
import numpy as np

def snptest(file):
    """ Read snptest data into an internal dataframe. """
    cols = ['chromosome','alleleA','alleleB','rsid','position','all_total', 'cases_total','controls_total','all_maf','frequentist_add_pvalue',
    'frequentist_add_beta_1', 'frequentist_add_se_1']
    df = pd.read_csv(file, sep=' ', comment='#')[cols]
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
    # takes MAF_U as reflects 'unaffected' population controls, unless MAF_U > 0.5, in which case it uses MAF_A
    frq_df.rename(columns={'CHR':'chromosome','SNP':'rsid','A2':'allele2','MAF_U':'maf','NCHROBS_A':'cases_total','NCHROBS_U':'controls_total'}, inplace=True)
    for index, row in frq_df.iterrows():
        if row['maf'] >= 0.5:
            row['maf'] = row['MAF_A']
    # if chromosome column has more than 1 number, read chromosome number from .assoc.logistic file and only include rows with that value.
    chromosomes = df.chromosome.unique()
    frq_df = frq_df[frq_df['chromosome'].isin(chromosomes)]
    frq_df = frq_df.drop(columns={"chromosome", "MAF_A"}, axis=1)
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
        map_file = pd.read_csv(file, sep='\t')
        chromosome = map_file['Chromosome'].ix[0].strip('chr')
        maps[chromosome] = map_file
    return maps

def annovar(file, file_exonic, colnames):
    """ Read ANNOVAR output files into an internal dataframe.

    Gene annotation with ANNOVAR returns two different output files (variant_function and exonic_variant_function).

    Where exonic SNPs exist, we merge the additional data of exonic variant function, and genes + transcript ID + protein-level change into the dataframe based on matching rsids.
    """
    df = pd.DataFrame(columns=colnames)
    if os.path.getsize(file) != 0:
        df = df.append(pd.read_csv(file, sep='\t', names = colnames))
    if os.path.getsize(file_exonic) != 0:
        df2 = pd.read_csv(file_exonic, sep='\t', names = colnames, usecols=range(1,(len(colnames) + 1)))
        df2 = df2.filter(items=['var_effect','genes','rsid'], axis=1)
        df2.rename(columns={'var_effect':'exonic_variant_function', 'genes':'genes_transcriptID'}, inplace=True)
        df2 = df2.set_index('rsid')
        df = pd.merge(df, df2, how='left',on='rsid')
    return df

def index(file):
    """ Read CRAFT .index output file into a dataframe."""
    index_df = pd.read_csv(file, sep='\t')
    return index_df

def abf_cred(file):
    """Read CRAFT .abf.cred file into a dataframe."""
    cred_snps = pd.read_csv(file, sep='\t')
    return cred_snps

def finemap_cred(file):
    """Read FINEMAP .cred file into a dataframe."""
    cred_snps = pd.read_csv(file, sep=' ')
    no_cols = int((len(cred_snps.columns) - 2)/ 2)
    cred_dfs = []
    for i in range(no_cols):
        cred_df = cred_snps[cred_snps.columns[2*i + 1: 2*i + 3]]
        cred_df = cred_df.rename(columns={cred_df.columns[0]:"rsid", cred_df.columns[1]:"pp"})
        cred_df = cred_df.dropna(axis=0)
        cred_dfs.append(cred_df)
    return cred_dfs

def cred_annotated(file):
    """Read CRAFT .cred.annotated file into a dataframe."""
    cred_df = pd.read_csv(file, sep='\t')
    return cred_df

def ld(file):
    """ Read CRAFT .ld output file into a numpy array."""
    ld_array = np.loadtxt(file)
    return ld_array

def variant_file(file):
    """ Read CRAFT variant_file rsids into a list."""
    variant_df = pd.read_csv(file, sep=' ')
    return variant_df

def snp(file):
    """Read FINEMAP .snp file into a dataframe."""
    snp_df = pd.read_csv(file, sep=' ')
    snp_df.rename(columns={'prob' : 'pp'}, inplace=True)
    return snp_df
