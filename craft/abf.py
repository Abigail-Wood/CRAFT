import numpy as np
import pandas as pd

from scipy.stats import norm

def calc_abf(pval, maf, n, n_controls, n_cases):
    """ Calculate Approximate Bayes Factor.

        (Wakefield, 2009, Genet Epidemiol.)
        Based on Chris Wallace work
    Args:
        pval (float): GWAS p-value
        maf (float): Minor allele freq
        n (int): Sample size
        n_controls: Number of controls
        n_cases: Number of cases
    Returns:
        ABF
    """

    # Assert/set types
    pval = float(pval)
    maf = float(maf)
    n = int(n)
    n_controls = int(n_controls)
    n_cases = int(n_cases)

    # Calculate Z-score
    z = np.absolute(norm.ppf(pval / 2))

    # Multiplicative model
    x0 = 0
    x1 = 1
    x2 = 2
    scale0 = n_controls
    scale1 = n_cases

    # No idea what this does
    d2 = (1-maf)**2 * x0**2 + 2*maf*(1-maf)*x1 + maf**2 * x2**2
    d1 = (1-maf)**2 * x0 + 2*maf*(1-maf)*x1 + maf**2 * x2
    V = (n_controls + n_cases) / (n_controls * n_cases * (d2-d1**2))

    # Scale V (by assignment of scaling parameters this is always 1)
    scale = ((n_controls + n_cases)/(scale0 + scale1)) * (scale0/n_controls) * (scale1/n_cases)
    V = V * scale

    # Compute W (no idea what this is either)
    W = (np.log(1.5) / norm.ppf(0.99))**2
    VW = V + W # simplification for ABF equation
    ABF = np.sqrt(VW/V) * np.exp(- z**2 * W / (2 * VW)) # 2 ln ABF, see Kass & Raftery 1995

    return ABF

def calc_postprob(data):
    """ Calculate posterior probability for each SNP."""
    sum_ABF = data['ABF'].sum()
    for index, row in data.iterrows():
        data['postprob'] = data['ABF'] / sum_ABF
    return data

def calc_postprobsum(data):
    """ Calc cumulative sum of the posterior probabilities."""
    for index, row in data.iterrows():
        data['postprob_cumsum'] = data['postprob'].cumsum()
    return data

def abf(data_dfs, cred_threshold):
    data_list = []
    for data in data_dfs:
        data['ABF'] = data.apply(
            lambda row: calc_abf(pval=row['pvalue'],
                                maf=row['maf'],
                                n=row['all_total'],
                                n_controls=row['controls_total'],
                                n_cases=row['cases_total']), axis=1)
        data = data.sort_values('ABF', ascending=False)
        data = calc_postprob(data)
        data = calc_postprobsum(data)
        data = data.sort_values('postprob_cumsum', ascending=False)
    # Trim credible SNPs based on posterior probability threshold
        if cred_threshold == '95':
            data = data[data.postprob_cumsum < 0.95]
        if cred_threshold =='99':
            data = data[data.postprob_cumsum < 0.99]
        data_list.append(data)
    return data_list
