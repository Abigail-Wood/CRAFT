import os
import re

import vcf as pyvcf

import craft.config as config

def base_annotation(df):
    """ TODO: summary docstring line.

    Annotates SNPs using VEP. Currently *not working*. Has dependencies.
    """
    # define file names
    annotation_file = "base_annotation.vep"
    annotation_vcf = "base_annotation.vcf"

    # create default VEP input
    annot_df = df.copy()
    annot_df['position2'] = annot_df['position']
    annot_df['alleles'] =  annot_df['alleleA'] + "/" + annot_df['alleleB']
    annot_df['strand'] = 1
    annot_df = annot_df[['chromosome','position','position2','alleles','strand']]
    annot_df.to_csv(annotation_file, sep='\t', index=False, header=False)

    # perform annotation with vep
    vep_cmd = "variant_effect_predictor.pl -i %s -o %s  --no_progress -q --assembly GRCh37 --dir_cache %s --symbol --force_overwrite --offline --pick --vcf --no_stats" % (annotation_file, annotation_vcf, config.ensembl_cache)
    os.system(vep_cmd)

    # format vcf
    vcf = pyvcf.Reader(open(annotation_vcf, 'r'))

    # define results table df
    info_fields = re.search('Format: (.+)', vcf.infos['CSQ'].desc).group(1).split("|")
    header = ['ID', 'CHROM', 'POS', 'REF', 'ALT'] + info_fields
    annot_results_df = pd.DataFrame(columns=header)

    for record in vcf:
        snp_info = [record.ID, record.CHROM, record.POS, record.REF, record.ALT]
        for csq in record.INFO['CSQ']:
            csq_info = csq.split("|")
            names = snp_info + csq_info
            result = pd.Series(names, index=header)
            annot_results_df = annot_results_df.append(result, ignore_index=True)

    # remove unwanted fields
    annot_results_df.CHROM = annot_results_df.CHROM.astype(int)
    annot_results_df.POS = annot_results_df.POS.astype(int)
    annot_results_df.rename(columns={'CHROM':'chromosome','POS':'position'}, inplace=True)
    annot_results_df = annot_results_df[['chromosome','position','Consequence','SYMBOL']]

    # merge annotation df with original df
    annot_results_df = pd.merge(df, annot_results_df)

    # clean up
    rm_cmd = "rm %s %s" % (annotation_file, annotation_vcf)
    os.system(rm_cmd)

    return annot_results_df

def merge_info(df, file):
    """ TODO: Docstring."""
    # read information file
    info = pd.read_table(file, delim_whitespace=True)

    # merge with df
    df = pd.merge(df, info, how='left', on='rsid')

    return df
