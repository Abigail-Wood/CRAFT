import sys
import os
import tempfile

import pandas as pd

import craft.config as config

def finemap(data_dfs, index_df, file_dir, n_causal_snps):
    """ Runs Finemap and LDStore on each SNP locus.

    Finemap(v1.3.1) was created by Christian Brenner (http://www.christianbenner.com/) and uses summary statistics for finemapping.

    **INPUT**
    1. Master file

    2. Z file (dataset, uncompressed)
    The dataset.z file is a space-delimited text file and contains the GWAS summary statistics one SNP per line. It contains the mandatory column names in the following order.

        `rsid` contains the SNP identifiers. The identifier can be a rsID number or a combination of chromosome name and genomic position (e.g. XXX:yyy)

        `chromosome` contains the chromosome names. The chromosome names can be chosen freely with precomputed SNP correlations (e.g. 'X', '0X' or 'chrX')

        `position` contains the base pair position for each SNP

        `allele1` contains the "first" allele of the SNPs. In SNPTEST this corresponds to 'allele_A', whereas BOLT-LMM uses 'ALLELE1'

        `allele2` contains the "second" allele of the SNPs. In SNPTEST this corresponds to 'allele_B', whereas BOLT-LMM uses 'ALLELE0'

        `maf` contains the minor allele frequencies

        `beta` contains the estimated effect sizes as given by GWAS software

        `se` contains the standard errors of effect sizes as given by GWAS software

    3. LD file
        Generated using LDstore, assuming PLINK input files (.bed, .bim, .fam).

    Other input options (support for BGEN input data, optional K file) are described at Christian Brenner's website and are not used here.

    **OUTPUT**


    """
    with tempfile.TemporaryDirectory() as tempdir:
        # make an empty master file
        master = pd.DataFrame(columns=['z','ld','snp','config','cred','log', 'n_samples'])

        ld_store_executable = os.path.join(config.ldstore_dir, "ldstore")
        master_file = os.path.join(tempdir, "master_file")
        master = open(master_file, "w")
        master.write("z;ld;snp;config;cred;log;n_samples\n")

        # need to take in index_df region definitions.
        index_count = 0

        for data in data_dfs:
            # set the PLINK basename based on chromosome in file
            chr = index_df.at[index_count, 'chromosome']
            index = index_df.at[index_count, 'rsid']

            # set filenames in tempdir with index SNP rsid (as unique identifier for input and output files)
            z_file = os.path.join(tempdir, index + ".z")
            variant_file = os.path.join(file_dir, index + "_variant.txt")
            plink_basename = os.path.join(config.plink_basename_dir, f"chr{chr}_ld_panel")
            bcor_file = os.path.join(tempdir, index + ".bcor")
            ld_file = os.path.join(file_dir, index + ".ld")
            snp_file = os.path.join(file_dir, index + ".snp")
            config_file = os.path.join(file_dir, index + ".config")
            cred_file = os.path.join(file_dir, index + ".cred")
            log_file = os.path.join(file_dir, index + ".log")

            # define region size [need to give index df as well and identify matching row based on rsid]
            region_start_cm = index_df.at[index_count, 'region_start_cm']
            region_end_cm = index_df.at[index_count, 'region_end_cm']

            # make a Z file
            order = ['rsid','chromosome','position','allele1','allele2','maf', 'beta', 'se']
            data = data[order]
            data.to_csv(z_file, sep=' ', index=False)

            # order of SNPs in LD file must correspond to order in Z file
            variants = data[['rsid','position','chromosome','allele1','allele2']]
            variants.to_csv(variant_file, sep=' ', index=False, header=['RSID','position','chromosome','A_allele','B_allele'])

            # make an LD file (bcor)
            cmd = (ld_store_executable + " --bplink " + plink_basename + f" --bcor {bcor_file} --incl-range {region_start_cm}-{region_end_cm} --n-threads 1")
            os.system(cmd)

            # make an LD file matrix for our rsids in locus (matrix)
            cmd = (ld_store_executable + f" --bcor {bcor_file}_1 --matrix {ld_file} --incl-variants {variant_file}")
            os.system(cmd)

            # append row to master file
            master.write(f"{z_file};{ld_file};{snp_file};{config_file};{cred_file};{log_file};{index_df.at[index_count, 'all_total']}\n")

            # increment index count to bring in new region definition.
            index_count+=1

        # Write completed master file out for use
        master.close()

        # run finemap (tell it data files are in temp directory)
        if n_causal_snps:
            cmd = (f"{config.finemap_dir}" + "/finemap_v1.3.1_x86_64" + f" --sss --in-files {master_file} --log  --n-causal-snps {n_causal_snps}")
            os.system(cmd)
        else:
            cmd = (f"{config.finemap_dir}" + "/finemap_v1.3.1_x86_64" + f" --sss --in-files {master_file} --log")
            os.system(cmd)
    return 0
