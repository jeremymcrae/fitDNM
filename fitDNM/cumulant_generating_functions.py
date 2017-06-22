
from math import isnan
from decimal import Decimal as dec

from numpy import exp

def cgf_0(mu, s, lambdas, weights):
    ''' cumulant generating function, see appendix A of Jiang et al 2015 AJHG 97:272
    
    This estimates the joint CGF of (Y,Z) as:
    
        \eqn{K_{Y,Z}(t,s) = \sum_{l=1}^{p} \lambda_l (e^{c_l t + s} - 1)}
    
    Args:
        mu: mutation rate
        s: something
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
    
    Returns:
        summed value
    '''
    
    assert len(lambdas) == len(weights)
    
    values = lambdas * (exp(weights * mu + s) - 1)
    
    total = sum(values)
    if not isnan(total):
        return total
    
    mu = dec(mu)
    s = dec(s)
    
    value = sum(( dec(m) * (exp(dec(w) * mu + s) - dec(1)) for m, w in zip(lambdas, weights) ))
    return float(value)

def cgf_2(mu, s, lambdas, weights):
    ''' cumulant generating function, second derivative
    
    The second derivative is calculated as:
    
    \eqn{K^{''}(t,s) = \begin{bmatrix}
          \sum_{l=1}^{p} \lambda_l {c_l}^2 (e^{c_l t + s} - 1) &
                 \sum_{l=1}^{p} \lambda_l {c_l} (e^{c_l t + s} - 1) \\
          \sum_{l=1}^{p} \lambda_l {c_l} (e^{c_l t + s} - 1) &
                 \sum_{l=1}^{p} \lambda_l (e^{c_l t + s} - 1) \\
        \end{bmatrix}}
    
    
    NOTE: The values calculated by this function differ very slightly from the
          original package function, due to floating point precision, given how
          the calculation is performed. This shouldn't change the overall result,
          and only makes a difference when sum(ss) * sum(mumu) is nearly equal to
          sum(mus)**2. Overall divergence should be less than 1e-17.
    
    Args:
        mu: mutation rate
        s: something
        lambdas: vector of per base and allele mutation rates
        weights: per base and allele weights, indicating the likely severity of
            each change
    
    Returns:
        determinant of the second derivative
    '''
    
    assert len(lambdas) == len(weights)
    
    standard = lambdas * exp(weights * mu + s)
    
    ss = standard
    mumu = weights**2 * standard
    mus = weights * standard
    
    if any(( isnan(x) for x in standard )):
        mu = dec(mu)
        s = dec(s)
        
        ss = dec(0)
        mumu = dec(0)
        mus = dec(0)
        for m, w in zip(lambdas, weights):
            m = dec(m)
            w = dec(w)
            
            x = m * exp(w - mu + s)
            ss += x
            mumu += w**2 * x
            mus += w * x
        
        total = ss * mumu - mus**2
        return float(total)
    else:
        return sum(ss) * sum(mumu) - (sum(mus)**2)
