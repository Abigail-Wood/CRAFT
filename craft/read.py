import glob
import pandas as pd

# Add cases_total or all_total for SNPtest - take from the column!

def snptest(file):
    """ Read snptest data into an internal dataframe. """
    cols = ['rsid','chromosome','position','alleleA','alleleB','info','cases_maf',
    'controls_maf','frequentist_add_pvalue']
    df = pd.read_table(file, sep=" ", comment="#")[cols]
    df.rename(columns={'frequentist_add_pvalue':'pvalue'}, inplace=True)
    return df

#
def plink(file, bim):
    """ Read PLINK (v1) data into an internal dataframe. """
    # read assoc. logistic
    cols = ['CHR','SNP','BP','P','OR','L95','U95']
    df = pd.read_table(file, delim_whitespace=True)[cols]
    df.rename(columns={'CHR':'chromosome','SNP':'rsid','BP':'position','A1':'alleleA','P':'pvalue'}, inplace=True)

    # read bim file
    cols = ['chromosome','rsid','cm','position','alleleA','alleleB']
    bdf = pd.read_table(bim, header=None, names=cols)
    bdf = bdf[['rsid','alleleA','alleleB']]

    # merge based on rsid
    df = pd.merge(df,bdf, how='inner',on='rsid')
    return df

def plink_noBIM(file):
    """ Read PLINK2 data into an internal  dataframe. """
    # read assoc.logistic
    cols = ['CHR','SNP','BP','A1','P', 'OR','L95','U95']
    df = pd.read_table(file, delim_whitespace=True)[cols]
    df.rename(columns={'CHR':'chromosome','SNP':'rsid','BP':'position','A1':'alleleA','P':'pvalue'}, inplace=True)
    df['alleleB'] = df['alleleA']
    return df

def indexsnps(file):
    """ Read CRAFT output table of index SNPs for fine-mapping use (TBD). """
    df = pd.read_table(file, sep="\t")
    return df

def generic(file):
    """Read tab-separated data into an internal dataframe. """
    df = pd.read_table(file, sep="\t")
    return df

def maps(source_dir):
    """ Read genetic map data into a maps object. """
    map_file_list = glob.glob(source_dir + "/*chr[0-9]*.txt")
    maps = {}
    for file in map_file_list:
        map_file = pd.read_table(file)
        chromosome = map_file['Chromosome'].ix[0].strip('chr')
        maps[chromosome] = map_file
    return maps
