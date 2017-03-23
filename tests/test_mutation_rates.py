# unit testing for the fitDNM functions

import unittest

import pandas

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from fitDNM.mutation_rates import get_gene_rates
from tests.compare_dataframes import CompareTables

class TestMutationRatesPy(CompareTables):
    ''' checks for obtaining table of mutation rates
    '''
    
    def test_get_gene_rates(self):
        ''' check that test_get_gene_rates is correct
        '''
        
        # pick a short gene to extra mutation rates for
        symbol = 'OR5A1'
        de_novos = pandas.DataFrame({'pos': [59211188], 'gene': [symbol]})
        
        mu_rates = get_gene_rates(symbol, de_novos)
        tab = StringIO('''gene chrom pos      ref alt           prob
                         OR5A1  11   59210642   A   C   1.931850e-09
                         OR5A1  11   59210642   A   G   1.199826e-08
                         OR5A1  11   59210642   A   T   2.196288e-09
                         OR5A1  11   59210643   T   A   2.130977e-09
                         OR5A1  11   59210643   T   C   1.203480e-08''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(mu_rates.head(), data)
