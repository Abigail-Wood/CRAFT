#!/usr/bin/env python

"""
Runs credible set analysis on a file of prefiltered summary stats.
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import norm

def parse_args():
    """ Load command line args."""
    parser = argparse.ArgumentParser()
    parser.add_argument('--p_col', metavar="<str>", help=('P-value column (default: pC)'), default="pC", type=str, required=False)
    parser.add_argument('--n_col', metavar="<str>", help=('N column (default: n)'), default="n", type=str, required=False)
    parser.add_argument('--maf_col', metavar="<str>", help=('MAF column (default: freq)'), default="freq", type=str, required=False)
    parser.add_argument('--prop_cases', metavar="<float>", help=('Proportion of cases, leave blank for quantitative traits (default: None)'), type=float, required=False)

    args = parser.parse_args()

    return args

def calc_abf(pval, maf, n, prop_cases=None):
    """ Calculate Approximate Bayes Factor (Wakefield, 2009, Genet Epidemiol.).
        Based on code from coloc: https://github.com/chr1swallace/coloc
    Args:
        pval (float): GWAS p-value
        maf (float): Minor allele freq
        n (int): Sample size
        n_cases (int or None): number of cases, if left blank will assume quantitative trait
    Returns:
        natural log(ABF)
    """

    # Chris Wallace based code - check comparisons

    # Assert/set types
    pval = float(pval)
    maf = float(maf)
    n = int(n)
    prop_cases = float(prop_cases) if prop_cases else None

    # Estimate variance for quant trait
    if prop_cases is None:
        sd_prior = 0.15
        v = var_data(maf, n)
    # Estimate var for cc study
    else:
        sd_prior = 0.2
        v = var_data_cc(maf, n, prop_cases)

    # Calculate Z-score
    z = np.absolute(norm.ppf(pval / 2))

    # Calc shrinkage factor: ratio of the prior variance to the total variance
    r = sd_prior**2 / (sd_prior**2 + v)

    # Approximate BF - ln scale to compare in log natural scale with LR diff
    lABF = 0.5 * (np.log(1 - r) + (r * z**2))

    return lABF

def log_sum(l):
    """ Calculates the log of the sum of the exponentiated logs taking out the
        max, i.e. insuring that the sum is not Inf
    Args:
        l (pandas Series)
    Returns:
        Sum of exponentiated logs
    """
    l_max = l.max()
    l_logsum = l_max + np.log(np.sum(np.exp(l - l_max)))
    return l_logsum

def freq_to_maf(freq):
    """ Convert allele frequency to minor allele freq
    """
    return min(freq, 1-freq)

def var_data(maf, n):
    """ Calc variance of MLE of beta for quantitative trait, assuming var(y)=0
    """
    var = 1 / (2 * n * maf * (1 - maf))
    return var

def var_data_cc(maf, n, prop_cases):
    """ Calc variance of MLE of beta for case-control
    """
    var = 1 / (2 * n * maf * (1 - maf) * prop_cases * (1 - prop_cases))
    return var
