# unit testing for the fitDNM functions

import unittest

from numpy import array

from fitDNM.cumulant_generating_functions import cgf_0, cgf_2

class TestConditionalApproximationPy(unittest.TestCase):
    ''' Cumulant generating function checks
    '''
    
    def test_cgf_0(self):
        ''' cgf_0 output is correct
        '''
        
        mu = 1
        s = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.05] * 6)
        self.assertAlmostEqual(cgf_0(mu, s, lambdas, weights), 1.11459067e-07, delta=1e-14)
        
        # expect an error if the vector lengths don't match
        lambdas = array([1e-8] * 5)
        with self.assertRaises(AssertionError):
            cgf_0(mu, s, lambdas, weights)
        
        # change the values and see if it still passes
        mu = 5
        s = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.5] * 6)
        self.assertAlmostEqual(cgf_0(mu, s, lambdas, weights), 1.926927118e-06, delta=1e-14)
        
        # now make a much longer vector and see if it still is correct
        mu = 1
        s = 1
        lambdas = array([1e-8] * 1000)
        weights = array([0.5] * 1000)
        self.assertAlmostEqual(cgf_0(mu, s, lambdas, weights), 3.4816890703380646611e-05, delta=1e-14)
    
    def test_cgf_2(self):
        ''' cgf_2 output is correct
        '''
        
        mu = 1
        s = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.05] * 6)
        self.assertAlmostEqual(cgf_2(mu, s, lambdas, weights), 1.232595164e-32, delta=1e-31)
        
        # expect an error if the vector lengths don't match
        lambdas = array([1e-8] * 5)
        with self.assertRaises(AssertionError):
            cgf_2(mu, s, lambdas, weights)
        
        # change the values and see if it still passes
        mu = 5
        s = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.5] * 6)
        self.assertEqual(cgf_2(mu, s, lambdas, weights), 0.0)
        
        # change the values and see if it still passes
        mu = 1
        s = 1
        lambdas = array([1e-4] * 100)
        weights = array([0.2] * 100)
        self.assertAlmostEqual(cgf_2(mu, s, lambdas, weights), 6.776264e-21, delta=1e-18)
