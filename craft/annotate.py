import os
import re
import tempfile

import vcf as pyvcf

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
        df.to_csv(to_annovar, sep='\t', index=False, header=False)
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

    FINEMAP outputs a .cred space-delimited text file. It contains the 95%
    credible sets for each causal signal conditional on other causal signals
    in the genomic region together with conditional posterior inclusion
    probabilities for each variant.

    ANNOVAR is able to take the list of rsids within the credible set and
    prepare it for input, using convert2annovar.pl

    ANNOVAR then adds gene-based annotation to the prepared input, using
    annotate_variation.pl.

    This function reads the final ANNOVAR input back and returns it as an
    internal dataframe, ready for merging with the original .cred file.
    """
    with tempfile.TemporaryDirectory() as tempdir:
        # make a list of rsids in credible SNP set
        rsid_list = list(cred_snps[cred_snps.columns[1]])
        # select locus DF information about rsids in credible SNP set
        locus_df = locus_df[locus_df['rsid'].isin(rsid_list)]
        cred_snps = prepare_df_annoVar(locus_df)
        # make file in tempdir
        to_annovar = os.path.join(tempdir, "to_annovar")
        cred_snps.to_csv(to_annovar, sep='\t', index=False, header=False)
        # perform annotation with ANNOVAR (give input, standard output)
        cmd = (f"{config.annovar_dir}/annotate_variation.pl -geneanno "
            "-dbtype refGene -buildver hg19 "
            f"{to_annovar} {config.annovar_dir}/humandb/")
        os.system(cmd)
        # read back in my temp output files as a dataframe with column names
        colnames = ['var_effect', 'genes'] + list(locus_df.columns)
        df = read.annovar(to_annovar + ".variant_function",
        to_annovar + ".exonic_variant_function", colnames)
    return df
