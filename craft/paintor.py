def paintor(data_dfs, index_df):
    """ Runs PAINTOR V3.0 on summary statistics.

    Usage information available at the PAINTOR wiki. https://github.com/gkichaev/PAINTOR_V3.0/wiki/2.-Input-Files-and-Formats

    The CRAFT pipeline does not implement visualisation with CANVIS (as this requires Python 2.7, which is near end-of-life.)
    """
    with tempfile.TemporaryDirectory() as tempdir:
        tempdir = "output/paintor_input/"

        ld_store_executable = os.path.join(config.ldstore_dir, "ldstore")
        input_file = os.path.join(tempdir, "input.file")
        input_file = open(input_file, "w")

        # need to take in index_df region definitions.
        index_count = 0

        for data in data_dfs:
            # set the PLINK basename based on chromosome in file
            chr = index_df.at[index_count, 'chromosome']
            index = index_df.at[index_count, 'rsid']

            # set filenames in tempdir with index SNP rsid (as unique identifier for input and output files)
            locus_file = os.path.join(tempdir, index)
            variant_file = os.path.join(tempdir, index + "_variant.txt")
            plink_basename = os.path.join(config.plink_basename_dir, f"chr{chr}_ld_panel")
            bcor_file = os.path.join(tempdir, index + ".bcor")
            ld_file = os.path.join(tempdir, index + ".ld")
            annotation_file = os.path.join(tempdir, index + ".annotations")

            # define region size [need to give index df as well and identify matching row based on rsid]
            region_start_cm = index_df.at[index_count, 'region_start_cm']
            region_end_cm = index_df.at[index_count, 'region_end_cm']

            # make and write a locus file
            order = ['chromosome','position','rsid', 'beta', 'se', 'allele1','allele2']
            data = data[order]
            # Z-score = beta / se (the Wald statistic)
            data['ZSCORE'] = data['beta']/data['se']
            data = data.drop(['beta','se'], axis=1)
            data.to_csv(locus_file, sep=' ', header=['CHR','POS','RSID','ZSCORE','ALLELE1','ALLELE2'])

            # order of SNPs in LD file must correspond to order in Z file
            variants = data[['rsid','position','chromosome','allele1','allele2']]
            variants.to_csv(variant_file, sep=' ', index=False, header=['RSID','position','chromosome','A_allele','B_allele'])

            # make an LD file (bcor)
            cmd = (ld_store_executable + " --bplink " + plink_basename + f" --bcor {bcor_file} --incl-range {region_start_cm}-{region_end_cm} --n-threads 1")
            os.system(cmd)

            # make an LD file matrix for our rsids in locus (matrix)
            cmd = (ld_store_executable + f" --bcor {bcor_file}_1 --matrix {ld_file} --incl-variants {variant_file}")
            os.system(cmd)

            # Make an annotation file (all rows 0 to show 'no annotation')
            # Annotation library (large, 6.7GB download) is available from PAINTOR and may be implemented in future versions of this pipeline
            annotation_df = pd.DataFrame(index=data.index,columns='dummy_annotation')
            annotation_df['dummy_annotation'] = 0
            annotation_df.to_csv(annotation_file, sep=' ', header=['dummy_annotation'])

            # append row to input file
            input_file.write(f"{index_df.at[index_count, 'rsid']}\n")

            # increment index count to bring in new region definition.
            index_count+=1

        # Write completed input file out for use
        input_file.close()

        # run paintor (tell it data files are in temp directory)
        # may wish to add command line option for specifying max causal and enumerate [number of causals]
        cmd = (f"{config.paintor_dir}" + "PAINTOR " + f" --input {input.file} -Zhead ZSCORE -LDname ld -in {tempdir} -out output/paintor_output - max_causal 2 -enumerate 2")

        os.system(cmd)

    return 0
