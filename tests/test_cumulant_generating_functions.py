# unit testing for the fitDNM functions

library(fitDNM)
library(testthat)

context("Cumulant generating function checks")

test_that("cgf_0 output is correct", {
    
    mu = 1
    s = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.05, 6)
    expect_equal(cgf_0(mu, s, lambdas, weights), 1.11459067e-07, tol=1e-14)
    
    # expect an error if the vector lengths don't match
    lambdas = rep(1e-8, 5)
    expect_error(cgf_0(mu, s, lambdas, weights))
    
    # change the values and see if it still passes
    mu = 5
    s = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.5, 6)
    expect_equal(cgf_0(mu, s, lambdas, weights), 1.926927118e-06, tol=1e-14)
    
    # now make a much longer vector and see if it still is correct
    mu = 1
    s = 1
    lambdas = rep(1e-8, 1000)
    weights = rep(0.5, 1000)
    expect_equal(cgf_0(mu, s, lambdas, weights), 3.4816890703380646611e-05, tol=1e-14)
})

test_that("cgf_2 output is correct", {
    
    mu = 1
    s = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.05, 6)
    expect_equal(cgf_2(mu, s, lambdas, weights), 1.232595164e-32, tol=1e-31)
    
    # expect an error if the vector lengths don't match
    lambdas = rep(1e-8, 5)
    expect_error(cgf_2(mu, s, lambdas, weights))
    
    # change the values and see if it still passes
    mu = 5
    s = 1
    lambdas = rep(1e-8, 6)
    weights = rep(0.5, 6)
    expect_identical(cgf_2(mu, s, lambdas, weights), 0)
    
    # change the values and see if it still passes
    mu = 1
    s = 1
    lambdas = rep(1e-4, 100)
    weights = rep(0.2, 100)
    expect_equal(cgf_2(mu, s, lambdas, weights), 6.776264e-21, tol=1e-18)
})
