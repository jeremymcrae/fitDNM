# unit testing for the fitDNM functions

import unittest

import pandas

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from fitDNM.mutation_rates import get_gene_rates, flatten_indexed_table
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
        tab = StringIO('''gene chrom pos      ref             A             C             G               T
                         OR5A1  11   59210642   A  0.000000e+00  1.931850e-09  1.199826e-08    2.196288e-09
                         OR5A1  11   59210643   T  2.130977e-09  1.203480e-08  1.920463e-09    0.000000e+00
                         OR5A1  11   59210644   G  9.635061e-09  2.477197e-09  0.000000e+00    2.657241e-09
                         OR5A1  11   59210645   T  1.816291e-09  5.247174e-09  1.499193e-09    0.000000e+00
                         OR5A1  11   59210646   C  2.085089e-09  0.000000e+00  2.086907e-09    7.784533e-09''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(mu_rates.head(), data)
    
    def test_flatten_indexed_table(self):
        '''
        '''
        
        # create a multi-index dataframe
        tuples = list(zip(*[['bar', 'bar', 'baz', 'baz',],
                        ['one', 'two', 'one', 'two',]]))
        
        index = pandas.MultiIndex.from_tuples(tuples, names=['first', 'second'])
        df = pandas.DataFrame({'A': range(0, 4), 'B': range(4, 8)}, index=index)
        
        # now create flattened version of the table
        tab = StringIO('''first second  A   B
                          bar   one     0   4
                          bar   two     1   5
                          baz   one     2   6
                          baz   two     3   7''')
        expected = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        expected = expected.reset_index(drop=True)
        
        self.compare_tables(flatten_indexed_table(df), expected)
