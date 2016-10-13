
from math import ceil

from scipy.stats import poisson

from fitDNM.conditional_approximation import conditional_approximation

def double_saddle_point_approximation(y, lambdas, weights):
    ''' estimate gene-wise de novo enrichment by saddlepoint approximation
    
    Args:
        y: this is the summed score for the observed de novos
        lambdas: vector of per base and allele mutation rates within a gene
        weights: vector of per base and allele weights, indicating the likely
            severity of each change.
    
    Returns:
        p-value for gene
    '''
    
    current = 0
    updated = 0
    total_mu = sum(lambdas)
    
    try:
        start = ceil(y / float(max(weights)))
    except ZeroDivisionError:
        return float('nan')
    
    if start <= 0:
        return 1
    
    # increment the expected score until the delta to the previous iteration is
    # less than one part in 10,000.
    for i in range(int(start), 101):
        updated = conditional_approximation(i, y, lambdas, weights) * poisson.pmf(i, total_mu)
        
        if updated is None:
            return None
        
        if abs(updated / current) < 1e-5:
            return current + updated
        
        current += updated
