import numpy as np
import pandas as pd

from scipy.stats import norm

def calc_abf(pval, maf, n, n_controls, n_cases):
    """Calculate Approximate Bayes Factor.

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

    **Usage**
    For a binary trait as no exact Bayes factor is calculable.
    For case-control studies only.

    **Limitations**
    Assumes causal variant is included in SNP population (i.e. requires relatively densely genotyped SNPs.)
    Assumes single causal variant only.
    """

    # Assert/set types
    pval = float(pval)
    maf = float(maf)
    n = int(n)
    n_controls = int(n_controls)
    n_cases = int(n_cases)

    # Variant of ABF calculation that uses p-values
    # Calculate Z-score
    z = np.absolute(norm.ppf(pval / 2))

    # Multiplicative model - assumes if 2 copies, 'twice as much' risk if have one copy.
    x0 = 0
    x1 = 1
    x2 = 2

    # No idea what this does
    d1 = (1-maf)**2 * x0 + 2*maf*(1-maf)*x1 + maf**2 * x2
    d2 = (1-maf)**2 * x0**2 + 2*maf*(1-maf)*x1 + maf**2 * x2**2
    V = (n_controls + n_cases) / (n_controls * n_cases * (d2-d1**2))

    # Compute W (no idea what this is either)
    # assumption: relative risk = 1.5
    # null hypothesis is that relative risk = 1
    # ppf = how many standard deviations is 99%
    W = (np.log(1.5) / norm.ppf(0.99))**2
    VW = V + W # simplification for ABF equation

    # Wakefield's approximate Bayes factor calculation (2009)
    ABF = np.sqrt(VW/V) * np.exp(- z**2 * W / (2 * VW))

    # Kass & Raftery (1995): 2 ln ABF allows comparison / rough interpretation of ABF meaning.
    # ln_2_ABF = 2 * np.log(np.sqrt(VW/V) * np.exp(- z**2 * W / (2 * VW)))

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
    #    if cred_threshold == '95':
    #        data = data[data.postprob_cumsum < 0.95]
    #    if cred_threshold =='99':
    #        data = data[data.postprob_cumsum < 0.99]
        data_list.append(data)
    return data_list
