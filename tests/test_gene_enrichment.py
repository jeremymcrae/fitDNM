# unit testing for the fitDNM functions

import unittest
import random

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

import pandas
from numpy import array
from numpy.random import beta, uniform, normal, seed

from fitDNM.gene_enrichment import compute_pvalue

class TestDoubleSaddlePointPy(unittest.TestCase):
    "check the primary gene enrichment function"
    
    def setUp(self):
        seed(1)
        random.seed(1)
    
    def test_compute_pvalue(self):
        ''' compute_pvalue output is correct
        '''
        
        tab = StringIO('''gene  chr  pos  ref  alt
                          GENE1   1    1    C    G
                          GENE1   1    2    A    T
                          GENEX   X    2    A    T''')
        
        n_male = 100
        n_female = 100
        de_novos = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        
        bases = ['A', 'C', 'G', 'T']
        severity = pandas.DataFrame({'gene': ['GENE1'] * 100,
            'chr': ['1'] * 100, 'pos': range(1, 101), 'ref': [ random.choice(bases) for x in range(100) ],
            'A': uniform(size=100), 'C': uniform(size=100),
            'G': uniform(size=100), 'T': uniform(size=100)})
        
        mu_rate = pandas.DataFrame({'gene': ['GENE1'] * 100,
            'chr': ['1'] * 100, 'pos': range(1, 101), 'ref': [ random.choice(bases) for x in range(100) ],
            'A': 10**(normal(size=100, loc=-8, scale=0.5)),
            'C': 10**(normal(size=100, loc=-8, scale=0.5)),
            'G': 10**(normal(size=100, loc=-8, scale=0.5)),
            'T': 10**(normal(size=100, loc=-8, scale=0.5))})
        
        symbol = 'GENE1'
        values = compute_pvalue(de_novos, n_male, n_female, symbol, severity, mu_rate)
        
        expected = {'symbol': 'GENE1', 'cohort_n': 200, 'nsnv_o': 196.0427494,
            'n_sites': 100, 'n_de_novos': 2, 'scores': 0.859,
            'p_value': 0.0003672926974, 'p_unweighted': 5.165670312e-06}
        
        self.assertEqual(values, expected)
