
from numpy import log, exp, sqrt
from scipy.stats import norm

from fitDNM.solver import solve_s_u
from fitDNM.cumulant_generating_functions import cgf_0, cgf_2

def get_w(s, s0, x, y, mu, lambdas, weights):
    ''' w part function
    
    Args:
        s: something
        s0: something
        x: tested summed score
        y: observed summed score
        mu: mutation rate
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
    
    Returns:
        summed value
    '''
    
    return 2 * (s * x + mu * y -
        cgf_0(mu=mu, s=s, lambdas=lambdas, weights=weights) -
        s0 * x + cgf_0(mu=0, s=s0 , lambdas=lambdas, weights=weights))

def avoid_zero_w_part(x, y, lambdas, weights, s0, increment=0.01):
    ''' adjust w_part until it is shifted away from zero
    
    Args:
        x: tested summed score
        y: observed summed score
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
        s0: initial s value
        increment: value to increment y value per iteration
    
    Returns:
        list of s, mu and w_part values
    '''
    while True:
        y += increment
        values = solve_s_u(x, y, lambdas, weights)
        values['w_part'] = get_w(values['s'], s0, x, y, values['mu'], lambdas, weights)
        
        if abs(values['w_part']) > 1e-4:
            break
    
    return values

def saddlepoint_p(values, s0, lambdas, weights):
    ''' calculate the p-value using the double saddle point
    
    Args:
        values: list of 's', 'mu' and 'w_part' values
        s0: initial s value
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
    
    Returns:
        p-value
    '''
    
    w = values['mu'] / abs(values['mu']) * sqrt(values['w_part'])
    
    # compute K_ss(s0, 0) and |K_2_(s, mu)|
    K_ss = exp(s0) * sum(lambdas)
    K_2 = cgf_2(values['mu'], values['s'], lambdas, weights)
    
    p_value = 1 - norm.cdf(w) + norm.pdf(w) * (sqrt(K_ss / (K_2)) / values['mu'] - 1 / w)
    
    return p_value

def approximate(x, y, lambdas, weights):
    ''' conditional approximation
    
    Args:
        x: tested summed score
        y: observed summed score
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
    
    Returns:
        approximate value
    '''
    
    # solve s0, s, mu and w_part
    s0 = log(x / sum(lambdas))
    values = solve_s_u(x, y, lambdas, weights)
    values['w_part'] = get_w(values['s'], s0, x, y, values['mu'], lambdas, weights)
    
    # ensure w_part >= 0
    i = 0
    while values['w_part'] < 0:
        # TODO: I haven't covered this section with a unit test yet
        i += 1
        values = solve_s_u(x, y, lambdas, weights, refine=10 * i)
        values['w_part'] = get_w(values['s'], s0, x, y, values['mu'], lambdas, weights)
        
        if i >= 10:
            return float('nan')
    
    if abs(values['w_part']) <= 1e-4:
        # If w_part is close enough to zero, this can throw off the estimate.
        # Avoid this by estimating w_part above and below, then use the average
        # NOTE: Possibly it only fails at exactly zero?
        lower = avoid_zero_w_part(x, y, lambdas, weights, s0, increment=0.01)
        lower_p = saddlepoint_p(lower, s0, lambdas, weights)
        
        upper = avoid_zero_w_part(x, y, lambdas, weights, s0, increment=-0.01)
        upper_p = saddlepoint_p(upper, s0, lambdas, weights)
        
        return (lower_p + upper_p) / 2
    else:
        return saddlepoint_p(values, s0, lambdas, weights)
