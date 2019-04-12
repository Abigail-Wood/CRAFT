import sys
import os

import pandas as pd
import numpy as np
import scipy.interpolate

def interpolate_cm(bp, map_file):
    """ TODO: summary docstring line.

    Given a base position and recombination map file return the recombination
    distance (in centimorgans, cM). Interpolate if required.
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

def interpolate_bp(cm, map_file):
    """ TODO: summary docstring line.

    Given a genetic distance (cM) and chromosome return the base position. Exact matching not performed as precision of floating points for cM a match is unlikely. Base position with closest cM is returned.
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

def get_index_snps_cm(df, alpha, distance, mhc, maps):
    """ Return a dataframe of index SNPs (with p > alpha)

    This function selects the SNP with the lowest p value > alpha, adds it to the list of index SNPs, discards everything within 'distance' range.
    """
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
        index_snp_cm = interpolate_cm(index_snp.position,maps[str(index_snp.chromosome.astype(int))])
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

def get_index_snps_bp(df, alpha, distance, mhc):
    """ TODO: docstring. """
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

def get_locus_snps(snps, index, distance_unit):
    """ Create dataframe of SNPs near index SNPs."""
    snps = snps.set_index('position')
    snps['index_rsid'] = ''
    locus_snps = pd.DataFrame(columns=snps.columns,index=pd.Index([],name='position'))
    data_dfs = []

    # For each index row, get the SNPs in that locus
    for index, row in index.iterrows():
        if distance_unit == 'cm':
            snps_in_locus = snps.loc[row.region_start_cm:row.region_end_cm].copy()
        elif distance_unit == 'bp':
            snps_in_locus = snps.loc[row.region_start_bp:row.region_end_bp].copy()
        snps_in_locus['index_rsid'] = row.rsid
        locus_snps = locus_snps.append(snps_in_locus)
        locus_snps = locus_snps.reset_index()
        # instead of appending snps_in_locus to locus_snps dataframe, add the dataframe to a list.
        data_dfs.append(locus_snps)
        # reset locus SNPs dataframe at the end of the loop.
        locus_snps = pd.DataFrame(columns=snps.columns,index=pd.Index([],name='position'))

    return data_dfs
