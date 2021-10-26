# unit testing for the fitDNM functions

import unittest
import random

from numpy import array
from numpy.random import beta, uniform, normal

from fitDNM.solver import solve_s_u

class TestSolverPy(unittest.TestCase):
    '''solver checks
    '''
    
    def test_solve_s_u(self):
        '''solve_s_u output is correct
        '''
        
        x = 1
        y = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.05] * 6)
        values = solve_s_u(x, y, lambdas, weights)
        
        self.assertAlmostEqual(values['mu'], 14195, delta=1e-14)
        self.assertAlmostEqual(values['s'], -693.12107872527565178, delta=1e-12)
        
        # check that the results are correct with different inputs.
        x = 2
        y = 1
        lambdas = array([1e-8] * 10)
        weights = array([0.1] * 10)
        values = solve_s_u(x, y, lambdas, weights)
        
        self.assertAlmostEqual(values['mu'], 7095, delta=1e-14)
        self.assertAlmostEqual(values['s'], -692.68875716848174307, delta=1e-12)
    
    def test_solve_s_u_with_delta(self):
        '''solve_s_u output is correct when we define a delta
        '''
        
        x = 1
        y = 1
        lambdas = array([1e-8] * 6)
        weights = array([0.05] * 6)
        values = solve_s_u(x, y, lambdas, weights, delta=200)
        
        self.assertAlmostEqual(values['mu'], 14100, delta=1e-14)
        self.assertAlmostEqual(values['s'], -688.37107872527565178, delta=1e-12)
    
    def test_solve_s_u_with_refining(self):
        '''solve_s_u output is correct when we need to refine
        '''
        
        x = 1
        y = 0.99
        lambdas = array([1e-6] * 100)
        weights = array([ (i + 1)/100.0 for i in range(100) ])
        start = 10
        values = solve_s_u(x, y, lambdas, weights, start=start)
        
        self.assertAlmostEqual(values['mu'], 69.31475000000003206, delta=1e-13)
        self.assertAlmostEqual(values['s'], -56.192386303155750227, delta=1e-13)
