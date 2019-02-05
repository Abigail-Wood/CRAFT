# module for craft_v2

import pandas as pd
import numpy as np
import scipy.interpolate
import sys
import os
import vcf as pyvcf
import re
import glob

import craft.config as config

def read_snptest(file):


	cols = ['rsid','chromosome','position','alleleA','alleleB','info','cases_maf','controls_maf','frequentist_add_pvalue']
	df = pd.read_table(file, sep=" ", comment="#")[cols]
	df.rename(columns={'frequentist_add_pvalue':'pvalue'}, inplace=True)

	return df

def read_generic(file):

	df = pd.read_table(file, sep="\t")

	return df

def read_plink(file, bim):

	# read assoc.logistic
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

def read_plink_2(file):

	# read assoc.logistic
	cols = ['CHR','SNP','BP','A1','P', 'OR','L95','U95']
	df = pd.read_table(file, delim_whitespace=True)[cols]
	df.rename(columns={'CHR':'chromosome','SNP':'rsid','BP':'position','A1':'alleleA','P':'pvalue'}, inplace=True)
	df['alleleB'] = df['alleleA']
	
	return df

def read_maps(source_dir):

	map_file_list = glob.glob(source_dir + "/*chr[0-9]*.txt")

	maps = {}

	for file in map_file_list:
		map_file = pd.read_table(file)
		chromosome = map_file['Chromosome'].ix[0].strip('chr')
		maps[chromosome] = map_file

	return maps

def merge_info(df, file):

	# read information file
	info = pd.read_table(file, delim_whitespace=True)

	# merge with df
	df = pd.merge(df, info, how='left', on='rsid')

	return df

def interpolate_cm(bp, map_file):

	"""

	Given a base position and recombination map file return the recombination distance (cM).
	Interpolate if required.

	"""

	if bp in map_file['Position(bp)'].values:
		cm = map_file['Map(cM)'].loc[map_file['Position(bp)'] == bp].real.astype(float)[0]

	else:
		index_location = np.searchsorted(np.array(map_file['Position(bp)']), bp)

		# get flanking base positions
		lower_bp = map_file['Position(bp)'].ix[index_location -1]
		upper_bp = map_file['Position(bp)'].ix[index_location]

		# interpolate map
		interpolation_range_bp = range(lower_bp,upper_bp,1)
		cm_interpolator = scipy.interpolate.interp1d(map_file['Position(bp)'], map_file['Map(cM)'])
		interpolation_range_cm = cm_interpolator(interpolation_range_bp)

		cm = interpolation_range_cm[interpolation_range_bp.index(bp)]

	return cm

def interpolate_bp(cm,map_file):

	"""

	Given a genetic distance (cM) and chromosome return the base position. Exact matching
	not performed as precision of floating points for cm a match is unlikely. Base position
	with closest cM is returned

	"""

	index_location = np.searchsorted(np.array(map_file['Map(cM)']), cm)

	# get flanking base positions
	lower_bp = map_file['Position(bp)'].ix[index_location -1]
	upper_bp = map_file['Position(bp)'].ix[index_location]

	# interpolate map
	interpolation_range_bp = range(lower_bp,upper_bp,1)
	cm_interpolator = scipy.interpolate.interp1d(map_file['Position(bp)'], map_file['Map(cM)'])
	interpolation_range_cm = cm_interpolator(interpolation_range_bp)

	bp = interpolation_range_bp[abs(interpolation_range_cm - cm).argmin()]

	return bp

def get_index_snps_bp(df, alpha, distance, mhc):

	# create df for results
	col_names = list(df.columns.values) + ['region_start_bp','region_end_bp']
	index_df = pd.DataFrame(columns=col_names)

	# exclude MHC region
	if not mhc:
		df = df.ix[~df.position.between(25000000, 35000000)]

	# get df of all index SNPs
	while df.pvalue.min() <= alpha:

		# get current index SNP
		index_snp = df.ix[df.pvalue.idxmin()]

		# define region boundaries
		region_start = index_snp.position - distance
		region_end = index_snp.position + distance
		index_snp = index_snp.append(pd.Series([region_start,region_end], index=['region_start_bp','region_end_bp']))

		# add current index SNP to index_df
		index_df = index_df.append(index_snp, ignore_index=True)

		# exclude index SNP region
		df = df.ix[~df.position.between(region_start, region_end)] 

	index_df.region_start_bp = index_df.region_start_bp.astype(int)
	index_df.region_end_bp = index_df.region_end_bp.astype(int)
	return index_df

def get_index_snps_cm(df, alpha, distance, mhc, maps):

	# create df for results
	col_names = list(df.columns.values) + ['region_start_cm','region_end_cm', 'region_size_kb']
	index_df = pd.DataFrame(columns=col_names)
	chr = df.chromosome.unique()[0]

	# exclude MHC region
	if not mhc and chr == 6:
		df = df.ix[~df.position.between(25000000, 35000000)]

	# get df of all index SNPs
	while df.pvalue.min() <= alpha:

		# get current index SNP
		index_snp = df.ix[df.pvalue.idxmin()]

		# define region boundaries
		index_snp_cm = interpolate_cm(index_snp.position, maps[str(index_snp.chromosome.astype(int))])
		region_start = interpolate_bp(index_snp_cm - distance,maps[str(index_snp.chromosome.astype(int))])
		region_end = interpolate_bp(index_snp_cm + distance,maps[str(index_snp.chromosome.astype(int))])
		region_size = round((region_end - region_start)/float(1000),1)
		index_snp = index_snp.append(pd.Series([region_start,region_end,region_size], index=['region_start_cm','region_end_cm','region_size_kb']))
		
		# add current index SNP to index_df
		index_df = index_df.append(index_snp, ignore_index=True)
		
		# exclude index SNP region
		df = df.ix[~df.position.between(region_start, region_end)] 

	index_df.region_start_cm = index_df.region_start_cm.astype(int)
	index_df.region_end_cm = index_df.region_end_cm.astype(int)
	return index_df

def base_annotation(df):

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
	#rm_cmd = "rm %s %s" % (annotation_file, annotation_vcf)
	#os.system(rm_cmd)

	return annot_results_df


