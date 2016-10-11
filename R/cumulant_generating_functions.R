
#' cumulant generating function, see appendix A of Jiang et al 2015 AJHG 97:272
#'
#' This estimates the joint CGF of (Y,Z) as:
#'
#'    \eqn{K_{Y,Z}(t,s) = \sum_{l=1}^{p} \lambda_l (e^{c_l t + s} - 1)}
#'
#' @param mu mutation rate
#' @param s something
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#'
#' @return summed value
cgf_0 <- function(mu, s, lambdas, weights) {
    
    stopifnot(length(lambdas) == length(weights))
    
    values = lambdas * (exp(weights * mu + s) - 1)
    
    return(sum(values))
}

#' cumulant generating function, second derivative
#'
#' The second derivative is calculated as:
#'
#' \eqn{K^{''}(t,s) = \begin{bmatrix}
#'      \sum_{l=1}^{p} \lambda_l {c_l}^2 (e^{c_l t + s} - 1) &
#'             \sum_{l=1}^{p} \lambda_l {c_l} (e^{c_l t + s} - 1) \\
#'      \sum_{l=1}^{p} \lambda_l {c_l} (e^{c_l t + s} - 1) &
#'             \sum_{l=1}^{p} \lambda_l (e^{c_l t + s} - 1) \\
#'    \end{bmatrix}}
#'
#'
#' NOTE: The values calculated by this function differ very slightly from the
#'       original package function, due to floating point precision, given how
#'       the calculation is performed. This shouldn't change the overall result,
#'       and only makes a difference when sum(ss) * sum(mumu) is nearly equal to
#'       sum(mus)**2. Overall divergence should be less than 1e-17.
#'
#' @param mu mutation rate
#' @param s something
#' @param lambdas vector of per base and allele mutation rates
#' @param weights per base and allele weights, indicating the likely severity of
#'        each change
#'
#' @return determinant of the second derivative
cgf_2 <- function(mu, s, lambdas, weights) {
    
    stopifnot(length(lambdas) == length(weights))
    
    standard = lambdas * exp(weights * mu + s)
    
    ss = standard
    mumu = weights**2 * standard
    mus = weights * standard
    
    return(sum(ss) * sum(mumu) - (sum(mus)**2))
}
