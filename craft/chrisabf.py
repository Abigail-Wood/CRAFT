import numpy as np
import pandas as pd

from scipy.stats import norm

def calc_abf(pval, maf, n, n_controls, n_cases):
    """ Calculate Approximate Bayes Factor (Wakefield, 2009, Genet Epidemiol.).
        Based on code from coloc: https://github.com/chr1swallace/coloc
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
    VW = V + W
    ABF = 2 * np.log(np.sqrt(VW/V)) ** (- z**2 * W / (2 * VW) )

    print(ABF)

    return ABF

def main():
    calc_abf(0.05, 0.01, 1000, 500, 500)

if __name__=='__main__':
    main()