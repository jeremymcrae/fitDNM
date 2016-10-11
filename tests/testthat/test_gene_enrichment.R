# unit testing for the fitDNM functions

library(fitDNM)
library(testthat)

context("check the primary gene enrichment function")

test_that("compute_pvalue output is correct", {
    
    set.seed(1)
    
    n_male = 100
    n_female = 100
    de_novos = read.table(text="
        gene  chr  pos  ref  alt
        GENE1   1    1    C    G
        GENE1   1    2    C    T
        GENEX   X    2    A    T", header=TRUE, stringsAsFactors=FALSE)
    
    bases = c('A', 'C', 'G', 'T')
    severity = data.frame('gene'=rep('GENE1', 100),
        'chr'=rep('1', 100), pos=seq(100), 'ref'=sample(bases, 100, replace=TRUE),
        'A'=runif(100), 'C'=runif(100), 'G'=runif(100), 'T'=runif(100))
    
    mu_rate = data.frame('gene'=rep('GENE1', 100),
        'chr'=rep('1', 100), pos=seq(100),'ref'=sample(bases, 100, replace=TRUE),
        'A'=10**(rnorm(100, mean=-8, sd=0.5)),
        'C'=10**(rnorm(100, mean=-8, sd=0.5)),
        'G'=10**(rnorm(100, mean=-8, sd=0.5)),
        'T'=10**(rnorm(100, mean=-8, sd=0.5)))
    
    symbol = 'GENE1'
    values = compute_pvalue(de_novos, n_male, n_female, symbol, severity, mu_rate)
    
    expected = list('symbol'='GENE1', 'cohort_n'=200, 'nsnv_o'=196.0427494,
        'n_sites'=100, 'n_de_novos'=2, 'scores'=0.859,
        'p_value'=0.0003672926974, 'p_unweighted'=5.165670312e-06)
    
    expect_equal(values, expected)
})
