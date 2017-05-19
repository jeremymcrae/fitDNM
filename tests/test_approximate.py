# unit testing for the fitDNM functions

import unittest
import random

from math import isnan

from numpy.random import beta, uniform, normal, seed

from fitDNM.approximate import approximate

class TestConditionalApproximationPy(unittest.TestCase):
    ''' unit test the conditional approximation function
    '''
    
    def setUp(self):
        ''' construct some objects for unit tests
        '''
        
        seed(1)
        random.seed(1)
    
    def test_approximate(self):
        ''' check approximate is correct
        '''
        
        x = 1
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = beta(a=0.5, b=0.5, size=1000)
        self.assertEqual(approximate(x, y, lambdas, weights), 0.011558509232964654)
    
    def test_approximate_low_w_part(self):
        ''' check approximate for low w_part values
        '''
        
        seed(13)
        
        x = 2
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = uniform(size=1000)
        self.assertEqual(approximate(x, y, lambdas, weights), 0.49792912429037084)
    
    def test_approximate_w_part_below_zero(self):
        ''' check approximate for w_part < 0
        '''
        
        x = 1.001
        y = 1
        lambdas = normal(loc=1e-4, scale=1e-5, size=1000)
        weights = uniform(size=1000)
        self.assertTrue(isnan(approximate(x, y, lambdas, weights)))
