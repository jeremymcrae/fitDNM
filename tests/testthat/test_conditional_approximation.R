# unit testing for the fitDNM functions

library(fitDNM)
library(testthat)

context("conditional approximation checks")

test_that("conditional_approximation output is correct", {
    
    set.seed(1)
    
    x = 1
    y = 1
    lambdas = rnorm(1000, mean=1e-4, sd=1e-5)
    weights = rbeta(1000, shape1=0.5, shape2=0.5)
    expect_equal(conditional_approximation(x, y, lambdas, weights), 0.01528782696)
})



test_that("conditional_approximation output is correct for low w_part values", {
    
    set.seed(1)
    
    x = 2
    y = 1
    lambdas = rnorm(1000, mean=1e-5, sd=1e-5)
    weights = runif(1000)
    expect_equal(conditional_approximation(x, y, lambdas, weights), 0.5095790724)
})
