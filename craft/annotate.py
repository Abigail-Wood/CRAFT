import os
import re
import tempfile

import vcf as pyvcf

import craft.config as config
import craft.read as read

def prepare_df_annoVar(df):
    """ TO-DO: docstring

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

def base_annotation_annoVar(df):
    """ TO-DO: docstring

    Docstring contents """
    with tempfile.TemporaryDirectory() as tempdir:
        # make a file in Temporary Directory, write to file
        to_annovar = os.path.join(tempdir, "to_annovar")
        df.to_csv(to_annovar, sep='\t', index=False, header=False)

        # perform annotation with ANNOVAR (give input, standard output)
        cmd = (f"{config.annovar_dir}/annotate_variation.pl -geneanno "
           "-dbtype refGene -buildver hg19 "
           f"{to_annovar} "
           f"{config.annovar_dir}/humandb/")
        os.system(cmd)
        # Output files written to -.variant_function, -.exonic_variant_function
        # add new columns and get original column names
        colnames = ['var_effect','genes'] + list(df.columns)
        # read back in my temp output files as a dataframe with column names
        df = read.annovar(to_annovar + ".variant_function", to_annovar + ".exonic_variant_function", colnames)
    return df
