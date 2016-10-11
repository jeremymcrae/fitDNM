# unit testing for the fitDNM functions

library(fitDNM)
library(testthat)

context("double saddle point approximation checks")

test_that("double_saddle_point_approximation output is correct", {
    
    set.seed(1)
    
    y = 1
    lambdas = rnorm(1000, mean=1e-4, sd=1e-5)
    weights = rbeta(1000, shape1=0.5, shape2=0.5)
    expect_equal(double_saddle_point_approximation(y, lambdas, weights), 0.002536151431)
})

test_that("double_saddle_point_approximation output is correct when the start
    value equals 0", {
    
    set.seed(1)
    
    # and check that if we have e
    y = 0
    lambdas = rnorm(1000, mean=1e-4, sd=1e-5)
    weights = rbeta(1000, shape1=0.5, shape2=0.5)
    expect_equal(double_saddle_point_approximation(y, lambdas, weights), 1)
    
    y = 1
    weights = rep(0, 1000)
    expect_equal(double_saddle_point_approximation(y, lambdas, weights), NA)
})



test_that("double_saddle_point_approximation output is correct", {
    
    set.seed(1)
    
    y = 1
    lambdas = rnorm(1000, mean=1e-5, sd=1e-5)
    weights = runif(1000)
    expect_equal(double_saddle_point_approximation(y, lambdas, weights), 2.477654823e-05)
})
