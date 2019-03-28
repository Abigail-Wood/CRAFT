import sys
import os
import craft.config as config

def finemap(data_dfs, ld_location):
    """ TO-DO: docstring

    Docstring contents """
    tempdir = "finemap_data"with tempfile.TemporaryDirectory() as tempdir:

        # make a z file
            # rsid chromosome position allele1 allele2 maf beta se
            # order of SNPS must correspond to order in LD
        for data in data_dfs
            Z_file = os.path.join(tempdir, "Z_file")
            df.to_csv(Z_file, sep='\t', index=False, header=False)
        # ensure LD file matches Z file SNP ordering
        # make a master file
            # z;ld;snp;config;cred;n_samples

        # perform annotation with ANNOVAR (give input, standard output)
        cmd = (f"{config.finemap_dir}" + "/finemap_v1.3.1_MacOSX")
        os.system(cmd)
        # Output files written to -.variant_function, -.exonic_variant_function
        # add new columns and get original column names
        colnames = ['var_effect','genes'] + list(df.columns)
        # read back in my temp output files as a dataframe with column names
        df = read.annovar(to_annovar + ".variant_function", to_annovar + ".exonic_variant_function", colnames)
return df
