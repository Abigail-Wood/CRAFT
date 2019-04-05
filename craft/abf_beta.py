import numpy as np
import pandas as pd

from scipy.stats import norm

def calc_abf_beta(beta, se, prior, log=False, log10=True):
    """ Wakefield ABF calculation

    Source: https://github.com/trochet/metabf

    beta is a single value or numeric vector. It should represent the observed effect size of a SNP on a trait from a genome-wide association study.

    se is a single value or numeric vector of the same length as beta. It should represent the standard error of the effect size of a SNP on a trait from a genome-wide association study.

    prior is a single value or numeric vector. If it is a single value, it does not need to be the same length as se and beta. This is the prior on true effect sizes.

    log sets whether the Bayes factor is returned in log (natural log) space.
    log10 sets whether the Bayes factor is returned in log10 space.
    """

    # Test of printing sum_abf and comparing with Finemap results indicates error is in ABF calculation

    #Setting up variables for calculation
    V = se**2 #var_beta
    W = prior**2 #SNPs are by default assumed to be causal with probability 1 / (# of SNPs in the genomic region)
    VW = V + W
    zsq = (beta**2) / V # Wald statistic

    if(log): #Perform calculation in natural log space
        log_abf = ((zsq/2)*(W/VW))-np.log(np.sqrt((VW)/V))
        return log_abf
    elif(log10): #Perform calculation in natural log space, convert to log10
        log10abf = np.log10(np.exp((zsq/2) * (W/VW))/np.sqrt(VW/V))
        return log10abf
    else:
        #Calculating Wakefield's approximate Bayes factor
        abf = np.exp((zsq/2) * (W/VW))/np.sqrt(VW/V)
        return abf

def calc_postprob(data):
    """ Calculate posterior probability for each SNP."""
    sum_ABF = data['log10ABF'].sum()
    print(sum_ABF)
    for index, row in data.iterrows():
        data['postprob'] = data['log10ABF'] / sum_ABF
    return data

def calc_postprobsum(data):
    """ Calc cumulative sum of the posterior probabilities."""
    for index, row in data.iterrows():
        data['postprob_cumsum'] = data['postprob'].cumsum()
    return data

def abf(data_dfs, cred_threshold):
    data_list = []
    for data in data_dfs:
        data['log10ABF'] = data.apply(
            lambda row: calc_abf_beta(beta=row['beta'],
                                se=row['se'],
                                prior=(1/len(data.index))), axis=1)
        data = data.sort_values('log10ABF', ascending=False)
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
