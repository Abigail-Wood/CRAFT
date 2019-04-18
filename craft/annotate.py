import os
import re
import tempfile

import vcf as pyvcf
import pandas as pd

import craft.config as config
import craft.read as read

def prepare_df_annoVar(df):
    """Prepare internal dataframe as input to ANNOVAR.

    Docstring contents """
    # make a list of all column names; position repeats twice for input
    df['position2'] = df['position']
    wanted = ['chromosome', 'position', 'position2','allele1', 'allele2']
    colnames = df.columns

    # list comprehensions to identify first 5 column names
    final_colnames = [col for col in wanted if col in colnames] + [col for col in colnames if col not in wanted]

    # re-order dataframe according to final list of column names and return
    annot_input = df[final_colnames]
    return annot_input

def annotation_annoVar(df):
    """Use ANNOVAR to annotate prepared internal dataframe.

    Describe ANNOVAR functions here. """
    with tempfile.TemporaryDirectory() as tempdir:
        # make file in tempdir, write to file
        to_annovar = os.path.join(tempdir, "to_annovar")
        df.to_csv(to_annovar, sep='\t', index=False, header=False, float_format='%g')
        # perform annotation with ANNOVAR (give input, standard output)
        cmd = (f"{config.annovar_dir}/annotate_variation.pl -geneanno "
           "-dbtype refGene -buildver hg19 "
           f"{to_annovar} {config.annovar_dir}/humandb/")
        os.system(cmd)
        # Output files written to -.variant_function, -.exonic_variant_function
        # add new columns and get original column names
        colnames = ['var_effect','genes'] + list(df.columns)
        # read back in my temp output files as a dataframe with column names
        df = read.annovar(to_annovar + ".variant_function",
        to_annovar + ".exonic_variant_function", colnames)
    return df

def finemap_annotation_annoVar(cred_snps, locus_df):
    """Use ANNOVAR to annotate prepared .cred FINEMAP output.

    FINEMAP outputs a .cred space-delimited text file. It contains the
    95% credible sets for each causal signal conditional on other causal
    signals in the genomic region together with conditional posterior
    inclusion probabilities for each variant.

    This function:
    Filters the locus dataframe (containing all summary statistic
    information, including chromosome, position, allele 1, allele2),
    using a list of rsids obtained from the .cred file.

    Uses the prepare_df_annoVAR function from craft.annotate to process
    the new dataframe as ANNOVAR input.

    Uses ANNOVAR to add gene-based annotation to the prepared input,
    using annotate_variation.pl.

    Reads the final ANNOVAR input back and returns it as an internal
    dataframe, removes unnecessary columns, and merges it using rsid as
    an index to add the posterior probability from the original  .cred
    file.
    """
    with tempfile.TemporaryDirectory() as tempdir:
        # make a list of rsids in credible SNP set
        rsid_list = list(cred_snps[cred_snps.columns[0]])
        # select locus DF information about rsids in credible SNP set
        locus_df = locus_df[locus_df['rsid'].isin(rsid_list)]
        cred_snps_prepared = prepare_df_annoVar(locus_df)
        # make file in tempdir
        to_annovar = os.path.join(tempdir, "to_annovar")
        cred_snps_prepared.to_csv(to_annovar, sep='\t', index=False,
                                  header=False, float_format='%g')
        # perform annotation with ANNOVAR (give input, standard output)
        cmd = (f"{config.annovar_dir}/annotate_variation.pl -geneanno "
               "-dbtype refGene -buildver hg19 "
               f"{to_annovar} {config.annovar_dir}/humandb/")
        os.system(cmd)
        # read back in my temp output files as a dataframe with column names
        colnames = ['var_effect', 'genes', 'chromosome','position','position2',
                    'allele1','allele2','rsid','all_total','cases_total',
                    'controls_total','maf','pvalue','beta','se','index_rsid',
                    'ABF','pp']
        df = read.annovar(to_annovar + ".variant_function",
        to_annovar + ".exonic_variant_function", colnames)
        # Drop unnecessary columns from locus SNPs dataframe before merge
        df = df.drop(['position2', 'ABF','pp'], axis=1)
        cred_snps = cred_snps.set_index('rsid')
        df = pd.merge(df, cred_snps, how='left',on='rsid')
        df = df.sort_values('pp', ascending=False)
    return df
