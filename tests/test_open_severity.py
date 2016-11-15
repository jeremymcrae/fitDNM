# unit testing for the fitDNM functions

import unittest
import os

import pandas

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from fitDNM.open_severity import get_cadd_severity
from tests.compare_dataframes import CompareTables

class TestOpenSeverityPy(CompareTables):
    ''' checks for obtaining table of severity scores
    '''
    
    def test_get_cadd_severity(self):
        ''' check that get_cadd_severity is correct
        '''
        
        # pick a short gene to extra mutation rates for
        symbol = 'OR5A1'
        chrom, start, end = '11', 59210641, 59211589
        
        cadd_path = os.path.join(os.path.dirname(__file__), 'data', 'cadd.txt.gz')
        
        severity = get_cadd_severity(symbol, chrom, start, end, cadd_path)
        tab = StringIO('''gene chrom pos      ref            A           C           G           T
                         OR5A1  11   59210642   A   0.00000000  0.86871956  0.53516373  0.87390438
                         OR5A1  11   59210643   T   0.89641425  0.10339738  0.84070578  0.00000000
                         OR5A1  11   59210644   G   0.89287273  0.84281901  0.00000000  0.86401218
                         OR5A1  11   59210645   T   0.00069053  0.00275928  0.00023023  0.00000000
                         OR5A1  11   59210646   C   0.08924828  0.00000000  0.82461194  0.35597904''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(severity.head(), data)
