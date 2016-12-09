
from __future__ import division

from math import isnan

from numpy import sign, exp, log

def cgf_ratio(lambdas, weights, mu):
    ''' get CGF-like ratio
    
    Args:
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
        mu: overall estimated mutation rate
    
    Returns:
        float value for
    '''
    return (sum(lambdas * weights * exp(weights * mu)) /
        sum(lambdas * exp(weights * mu)))

def solver(current, initial_sign, ratio, lambdas, weights, delta):
    ''' solve for initial and updated
    
    Args:
        current: currently estimated mutation rate
        initial_sign: sign of initial estimate
        ratio: ratio of estimated and observed summed scores
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
        delta: amount to update by during iterations
    
    Returns:
        dictionary of inital and updated values
    '''
    refine = False
    while True:
        updated = current + initial_sign * delta
        value = ratio - cgf_ratio(lambdas, weights, updated)
        
        if isnan(sign(value)):
            # TODO: check that this clause is correct. Perhaps it should check
            # TODO: if value is None, or sign(value) == -1L?
            break
        elif initial_sign != sign(value):
            refine = True
            break
        else:
            current = updated
    
    return {'current': current, 'updated': updated, 'refine': refine}

def solve_s_u(x, y, lambdas, weights, delta=10, refine=5, start=0):
    ''' solve for S and mu
    
    Args:
        x: estimated summed score
        y: observed summed score
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
        delta: amount to update by during iterations
        refine: number of iteratative refinements to perform
        start: initial mu estimate
    
    Returns:
        dictionary with 'mu' and 's' entries
    '''
    # solve mu first
    ratio = y / x
    initial_sign = sign(ratio - cgf_ratio(lambdas, weights, start))
    
    val = solver(start, initial_sign, ratio, lambdas, weights, delta)
    
    if val['refine']:
        for i in range(refine):
            delta = delta / 10
            val = solver(val['current'], initial_sign, ratio, lambdas, weights, delta)
    
    mu = (val['current'] + val['updated']) / 2
    s = log(x) - log(sum(lambdas * exp(weights * mu)))
    
    return {'s': s, 'mu': mu}
