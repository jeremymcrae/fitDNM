
#' get CGF-like ratio
#'
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#' @param mu overall estimated mutation rate
#'
#' @return float value for
cgf_ratio <- function(lambdas, weights, mu) {
    return(sum(lambdas * weights * exp(weights * mu)) /
        sum(lambdas * exp(weights * mu)))
}

#' solve for initial and updated
#'
#' @param current currently estimated mutation rate
#' @param initial_sign sign of initial estimate
#' @param ratio ratio of estimated and observed summed scores
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#' @param delta amount to update by during iterations
#'
#' @return list of inital and updated values
solver <- function(current, initial_sign, ratio, lambdas, weights, delta) {
    refine = FALSE
    while (TRUE) {
        updated = current + initial_sign * delta
        value = ratio - cgf_ratio(lambdas, weights, updated)
        
        if (is.na(sign(value))) {
            break
        } else if (initial_sign != sign(value)) {
            refine = TRUE
            break
        } else {
            current = updated
        }
    }
    
    return(list('current'=current, 'updated'=updated, 'refine'=refine))
}

#' solve for S and mu
#'
#' @param x estimated summed score
#' @param y observed summed score
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#' @param delta amount to update by during iterations
#' @param refine number of iteratative refinements to perform
#' @param start initial mu estimate
#'
#' @return list with 'mu' and 's' entries
solve_s_u <- function(x, y, lambdas, weights, delta=10, refine=5, start=0) {
    # solve mu first
    ratio = y / x
    initial_sign = sign(ratio - cgf_ratio(lambdas, weights, start))
    
    val = solver(start, initial_sign, ratio, lambdas, weights, delta)
    
    if (val$refine) {
        for (i in seq(refine)) {
            delta = delta / 10
            val = solver(val$current, initial_sign, ratio, lambdas, weights, delta)
        }
    }
    
    mu = (val$current + val$updated) / 2
    s = log(x / sum(lambdas * exp(weights * mu)))
    
    return(list('s'=s, 'mu'=mu))
}
