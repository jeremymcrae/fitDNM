# unit testing for the fitDNM functions

import unittest
import random
from math import isnan

from numpy import array
from numpy.random import beta, uniform, normal, seed

from fitDNM.saddlepoint import saddlepoint

class TestDoubleSaddlePointPy(unittest.TestCase):
    ''' double saddle point approximation checks
    '''
    
    def setUp(self):
        seed(1)
        random.seed(1)
    
    def test_saddlepoint(self):
        '''saddlepoint output is correct
        '''
        
        y = 1
        lambdas = normal(size=1000, loc=1e-4, scale=1e-5)
        weights = beta(size=1000, a=0.5, b=0.5)
        self.assertAlmostEqual(saddlepoint(
            y, lambdas, weights), 0.0022951333364715112, delta=1e-12)
    
    def test_saddlepoint_zero_values(self):
        '''saddlepoint output is correct when the start value equals 0
        '''
        
        # and check that if we have e
        y = 0
        lambdas = normal(size=1000, loc=1e-4, scale=1e-5)
        weights = beta(size=1000, a=0.5, b=0.5)
        self.assertEqual(saddlepoint(y, lambdas, weights), 1)
        
        y = 1
        weights = array([0] * 1000)
        self.assertTrue(isnan(saddlepoint(y, lambdas, weights)))
    
    def test_saddlepoint_uniform(self):
        '''saddlepoint output is correct
        '''
        
        y = 1
        lambdas = normal(size=1000, loc=1e-5, scale=1e-5)
        weights = uniform(size=1000)
        self.assertEqual(saddlepoint(y, lambdas, weights), 2.7451994001564722e-05)
