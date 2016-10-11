
#' w part function
#'
#' @param s something
#' @param s0 something
#' @param x something
#' @param y something
#' @param mu mutation rate
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#'
#' @return summed value
get_w <- function(s, s0, x, y, mu, lambdas, weights) {
    
    return(2 * (s * x + mu * y -
        cgf_0(mu=mu, s=s, lambdas=lambdas, weights=weights) -
        s0 * x + cgf_0(mu=0, s=s0 , lambdas=lambdas, weights=weights)))
}

#' refine w part until abs(w_part) is > 1e-4
#'
#' @param x something
#' @param y something
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#' @param s0 initial s value
#' @param increment value to increment y value per iteration
#'
#' @return list of s, mu and w_part values
refine_w <- function(x, y, lambdas, weights, s0, increment=0.01) {
    while (TRUE) {
        y = y + increment
        values = solve_s_u(x, y, lambdas, weights)
        values['w_part'] = get_w(values$s, s0, x, y, values$mu, lambdas, weights)
        
        if (abs(values$w_part) > 1e-4) { break }
    }
    
    return(values)
}

#' calculate the p-value using the double saddle point
#'
#' @param values list of 's', 'mu' and 'w_part' values
#' @param s0 initial s value
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#'
#' @return p-value
saddlepoint_p <- function(values, s0, lambdas, weights) {
    
    w = values$mu / abs(values$mu) * sqrt(values$w_part)
    
    # compute K_ss(s0, 0) and |K_2_(s, mu)|
    K_ss_s0 = exp(s0) * sum(lambdas)
    K_2_smu = cgf_2(mu=values$mu, s=values$s, lambdas, weights)
    
    p_value = 1 - pnorm(w) + dnorm(w) * (sqrt(K_ss_s0 / (K_2_smu)) / values$mu - 1 / w)
    
    return(p_value)
}

#' conditional approximation
#'
#' @param x estimated summed score
#' @param y observed summed score
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#'
#' @return approximate value
conditional_approximation <- function(x, y, lambdas, weights) {
    
    # solve s0, s, mu and w_part
    s0 = log(x / sum(lambdas))
    values = solve_s_u(x, y, lambdas, weights)
    values['w_part'] = get_w(values$s, s0, x, y, values$mu, lambdas, weights)
    
    # ensure w_part >= 0
    refine = 10
    iter = 0
    while (values$w_part < 0 & iter < 10) {
        # TODO: I haven't covered this section with a unit test yet
        iter = iter + 1
        values = solve_s_u(x, y, lambdas, weights, refine=refine * iter)
        values['w_part'] = get_w(values$s, s0, x, y, values$mu, lambdas, weights)
    }
    
    if (iter >= 100) { return(NA) }
    
    if (abs(values$w_part) <= 1e-4) {
        lower = refine_w(x, y, lambdas, weights, s0, 0.01)
        lower_p = saddlepoint_p(lower, s0, lambdas, weights)
        
        upper = refine_w(x, y, lambdas, weights, s0, -0.01)
        upper_p = saddlepoint_p(upper, s0, lambdas, weights)
        
        return((lower_p + upper_p) / 2)
    } else {
        return(saddlepoint_p(values, s0, lambdas, weights))
    }
}
