# unit testing for the fitDNM functions

import unittest
import random

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
        self.assertEqual(conditional_approximation(x, y, lambdas, weights), 0.01528782696)
    
    def test_conditional_approximation_low_w_part(self):
        ''' check conditional_approximation for low w_part values
        '''
        
        x = 2
        y = 1
        lambdas = normal(loc=1e-5, scale=1e-5, size=1000)
        weights = uniform(size=1000)
        self.assertEqual(conditional_approximation(x, y, lambdas, weights), 0.01528782696)
