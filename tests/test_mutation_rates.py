# unit testing for the fitDNM functions

import unittest
import tempfile
import gzip

import os

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
        
        data = [['transcript', 'gene', 'chr', 'amino_acids', 'obs_exp', 'chisq_diff_null'],
            ['ENST00000302030.2', 'OR5A1', '1', '1-1', '0.1', '51'],
            ['ENST00000302030.2', 'OR5A1', '1', '2-315', '0.810418512627', '15.5']]
        data = '\n'.join([ '\t'.join(x) for x in data ])
        constraint = tempfile.NamedTemporaryFile(suffix='.txt.gz')
        
        with gzip.open(constraint.name, 'wt') as handle:
            handle.write(data)
        
        # pick a short gene to extra mutation rates for
        symbol = 'OR5A1'
        de_novos = pandas.DataFrame({'pos': [59211188], 'gene': [symbol]})
        
        mu_rates = get_gene_rates(symbol, de_novos, constraint.name)
        tab = StringIO('''gene chrom pos      ref alt           prob  consequence  constrained
                         OR5A1  11   59210642   A   C   1.931850e-09  missense     True
                         OR5A1  11   59210642   A   G   1.199826e-08  missense     True
                         OR5A1  11   59210642   A   T   2.196288e-09  missense     True
                         OR5A1  11   59210643   T   A   2.130977e-09  missense     True
                         OR5A1  11   59210643   T   C   1.203480e-08  missense     True
                         OR5A1  11   59210643   T   G   1.920463e-09  missense     True
                         OR5A1  11   59210644   G   A   9.635061e-09  missense     True
                         OR5A1  11   59210644   G   C   2.477197e-09  missense     True
                         OR5A1  11   59210644   G   T   2.657241e-09  missense     True
                         OR5A1  11   59210645   T   A   1.816291e-09  missense     False''')
        
        data = pandas.read_table(tab, sep='\s+', skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(mu_rates.head(10), data)
        
        mu_rates = get_gene_rates(symbol, de_novos)
        tab = StringIO('''gene chrom pos      ref alt           prob  consequence  constrained
                         OR5A1  11   59210642   A   C   1.931850e-09  missense     False
                         OR5A1  11   59210642   A   G   1.199826e-08  missense     False
                         OR5A1  11   59210642   A   T   2.196288e-09  missense     False
                         OR5A1  11   59210643   T   A   2.130977e-09  missense     False
                         OR5A1  11   59210643   T   C   1.203480e-08  missense     False
                         OR5A1  11   59210643   T   G   1.920463e-09  missense     False
                         OR5A1  11   59210644   G   A   9.635061e-09  missense     False
                         OR5A1  11   59210644   G   C   2.477197e-09  missense     False
                         OR5A1  11   59210644   G   T   2.657241e-09  missense     False
                         OR5A1  11   59210645   T   A   1.816291e-09  missense     False''')
        
        data = pandas.read_table(tab, sep='\s+', skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(mu_rates.head(10), data)
