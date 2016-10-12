#' exclude certain values from a table ,and convert to a vector
#'
#' @param data table of values
#' @param exclude matching table of booleans, where TRUE indicates values to exclude
#'
#' @export
#'
#' @return vector of values, with NA values removed
exclude_and_get_vector <- function(data, exclude) {
    data[exclude] = NA
    data = as.vector(data)
    
    return(data[!is.na(data)])
}

#' compute de novo enrichment for a gene
#'
#' @param de_novos table of de novos
#' @param n_male number of males
#' @param n_female number of females
#' @param symbol HGNC symbol for a gene
#' @param severity table of per base and per allele severity scores, for every gene
#' @param mu_rate table of mutation rates
#'
#' @export
#'
#' @return vector of p-value of de novo enrichment for a gene
compute_pvalue <- function(de_novos, n_male, n_female, symbol, severity, mu_rate) {
    
    # select the mutations for the given gene
    de_novos = de_novos[de_novos$gene == symbol, ]
    
    # count the samples in the cohort, so we can calculate expected mutations
    nsample = n_male + n_female
    if (de_novos$chr[1] == "X") {
        nsample = (n_male + n_female * 2) / 2
    }
    
    bases = c('A', 'C', 'G', 'T')
    
    severity = severity[severity$gene == symbol, ]
    positions = severity$pos
    ref_allele = severity$ref
    severity = severity[, bases]
    nsnv_o = sum(severity)
    mu_rate = mu_rate[mu_rate$gene == symbol, bases]
    
    # recode for LOF mutations
    severity[severity == 2] = 1
    severity[severity == 3] = 0
    
    if (nrow(severity) == 0) { return(NULL) }
    
    # downweight the severities by the rate across the alleles at each site
    severity = severity - apply((severity * mu_rate), 1, sum)
    
    # get the expected mutation rates per base per site for the gene
    mu_rate = nsample * 2 * mu_rate
    
    if (any(is.na(severity))) { warning(symbol, " has NA in severity") }
    
    # check if the position is in the CCDS (not included promoter region)
    de_novos = de_novos[de_novos$pos %in% positions, ]
    
    nsnv = nrow(severity)
    ndenovo = nrow(de_novos)
    scores = 0
    
    # if there aren't any sites with mutation rates, or we don't have any de
    # novos, give the gene a p-value of 1
    if (nrow(mu_rate) <= 0 | ndenovo == 0) {
        p_value = 1
        p_unweighted = 1
    } else {
        for (idx in seq(nrow(de_novos))) {
            if (de_novos$ref[idx] != ref_allele[positions == de_novos$pos[idx]]) {
                stop(de_novos$chr[idx], de_novos$pos[idx], de_novos$ref[idx], "do not match ref\n")
            }
            scores = scores + severity[positions == de_novos$pos[idx], de_novos$alt[idx]]
        }
        
        # exclude data from non-mutated, i.e. unchanged sites
        exclude = mu_rate == 0
        mu_rate = exclude_and_get_vector(mu_rate, exclude)
        severity = exclude_and_get_vector(severity, exclude)
        
        # compute p-values, using saddlepoint and standard poisson approaches
        p_value = double_saddle_point_approximation(scores, mu_rate, severity)
        p_unweighted = 1 - ppois(nrow(de_novos) - 1, sum(mu_rate))
    }
    
    values = list('symbol'=symbol, 'cohort_n'=nsample, 'nsnv_o'=nsnv_o,
        'n_sites'=nsnv, 'n_de_novos'=ndenovo, 'scores'=round(scores, 3),
        'p_value'=p_value, 'p_unweighted'=p_unweighted)
    
    return(values)
}
