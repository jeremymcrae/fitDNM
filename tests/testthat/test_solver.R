# unit testing for the fitDNM functions

library(fitDNM)
library(testthat)

context("solver checks")

test_that("solve_s_u output is correct", {
    
    x = 1
    y = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.05, 6)
    values = solve_s_u(x, y, lambdas, weights)
    
    expect_equal(values$mu, 14195, tol=1e-14)
    expect_equal(values$s, -693.12107872527565178, tol=1e-14)
    
    # check that the results are correct with different inputs.
    x = 2
    y = 1
    lambdas = rep(1e-8, 10)
    weights = rep(0.1, 10)
    values = solve_s_u(x, y, lambdas, weights)
    
    expect_equal(values$mu, 7095, tol=1e-14)
    expect_equal(values$s, -692.68875716848174307, tol=1e-14)
})

test_that("solve_s_u output is correct when we define a delta", {
    
    x = 1
    y = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.05, 6)
    values = solve_s_u(x, y, lambdas, weights, delta=200)
    
    expect_equal(values$mu, 14100, tol=1e-14)
    expect_equal(values$s, -688.37107872527565178, tol=1e-14)
})


test_that("solve_s_u output is when we need to refine", {
    
    x = 1
    y = 0.99
    lambdas = rep(1e-6, 100)
    weights = seq(100)/100
    start = 10
    values = solve_s_u(x, y, lambdas, weights, start=start)
    
    expect_equal(values$mu, 69.31475000000003206, tol=1e-14)
    expect_equal(values$s, -56.192386303155750227, tol=1e-14)
})
