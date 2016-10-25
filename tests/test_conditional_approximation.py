# unit testing for the fitDNM functions

import unittest
import random

from math import isnan

from numpy.random import beta, uniform, normal, seed

from fitDNM.conditional_approximation import conditional_approximation

class TestConditionalApproximationPy(unittest.TestCase):
    ''' unit test the conditional approximation function
    '''
    
    def setUp(self):
        ''' construct some objects for unit tests
        '''
        
        seed(1)
        random.seed(1)
    
    def test_conditional_approximation(self):
        ''' check conditional_approximation is correct
        '''
        
        x = 1
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = beta(a=0.5, b=0.5, size=1000)
        self.assertEqual(conditional_approximation(x, y, lambdas, weights), 0.011558509232964654)
    
    def test_conditional_approximation_low_w_part(self):
        ''' check conditional_approximation for low w_part values
        '''
        
        seed(13)
        
        x = 2
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = uniform(size=1000)
        self.assertEqual(conditional_approximation(x, y, lambdas, weights), 0.49792912429037084)
    
    def test_conditional_approximation_w_part_below_zero(self):
        ''' check conditional_approximation for w_part < 0
        '''
        
        x = 1.001
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = uniform(size=1000)
        self.assertTrue(isnan(conditional_approximation(x, y, lambdas, weights)))