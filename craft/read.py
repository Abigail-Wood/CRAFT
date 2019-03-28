import glob
import os

import pandas as pd

def snptest(file):
    """ Read snptest data into an internal dataframe. """
    cols = ['chromosome','alleleA','alleleB','rsid','position','all_total', 'cases_total','controls_total','all_maf','frequentist_add_pvalue',
    'frequentist_add_beta_1', 'frequentist_add_se_1']
    df = pd.read_table(file, sep=' ', comment='#')[cols]
    df.rename(columns={'frequentist_add_pvalue':'pvalue', 'frequentist_add_beta_1':'beta', 'frequentist_add_se_1':'se'}, inplace=True)
    return df

def plink(file, bim):
    """ Read PLINK (v1) data into an internal dataframe. """
    # read assoc. logistic
    cols = ['CHR','SNP','BP','SE','P','OR','L95','U95']
    df = pd.read_table(file, delim_whitespace=True)[cols]
    df.rename(columns={'CHR':'chromosome','SNP':'rsid','BP':'position','A1':'alleleA','P':'pvalue', 'SE':'se'}, inplace=True)

    # Add column for all cases, all total, all control, minor allele frequency, beta. Read from BIM file?

    # read bim file
    cols = ['chromosome','rsid','cm','position','alleleA','alleleB']
    bdf = pd.read_table(bim, header=None, names=cols)
    bdf = bdf[['rsid','alleleA','alleleB']]

    # merge based on rsid
    df = pd.merge(df,bdf, how='inner',on='rsid')
    return df

def indexsnps(file):
    """ Read in CRAFT output file of index SNPs."""
    cols = ['chromosome','rsid','alleleA','alleleB','position','all_total', 'cases_total','controls_total','all_maf','pvalue',
    'beta', 'se','region_start_cm','region_end_cm','region_size_kb']
    df = pd.read_table(file, sep='\t')
    return df

def csv(file):
    """Read csv data into an internal dataframe. """
    cols = ['chromosome','rsid','alleleA','alleleB','position','all_total', 'cases_total','controls_total','all_maf','pvalue',
    'beta', 'se']
    df = pd.read_table(file, sep='\t')
    return df


    df = pd.read_table(file, sep=' ', comment='#')[cols]
    df.rename(columns={'frequentist_add_pvalue':'pvalue', 'frequentist_add_beta_1':'beta', 'frequentist_add_se_1':'se'}, inplace=True)
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
