
#' estimate gene-wise de novo enrichment by saddlepoint approximation
#'
#' @param y this is the summed score for the observed de novos
#' @param lambdas vector of per base and allele mutation rates within a gene
#' @param weights vector of per base and allele weights, indicating the likely
#'        severity of each change.
#'
#' @export
#'
#' @return p-value for gene
double_saddle_point_approximation <- function(y, lambdas, weights) {
    
    current = 0
    updated = 0
    total_mu = sum(lambdas)
    
    max_weights = max(weights, na.rm=TRUE)
    start = ceiling(y / max_weights)
    
    if (start <= 0) {
        return(1)
    } else if (start == Inf) {
        return(NA)
    }
    
    # increment the expected score until the delta to the previous iteration is
    # less than one part in 10,000.
    for (i in c(start:100)) {
        updated = conditional_approximation(i, y, lambdas, weights) * dpois(i, total_mu)
        
        if (is.na(updated)) { return(NA) }
        
        if (abs(updated / current) < 1e-5) {
            return(current + updated)
        }
        
        current = current + updated
    }
}
