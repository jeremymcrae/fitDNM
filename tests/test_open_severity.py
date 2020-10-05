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
        tab = StringIO('''gene chrom pos      ref  alt   score
                         OR5A1  11   59210642   A    C   8.818
                         OR5A1  11   59210642   A    G   3.327
                         OR5A1  11   59210642   A    T   8.993
                         OR5A1  11   59210643   T    A   9.847
                         OR5A1  11   59210643   T    C   0.474''')
        
        data = pandas.read_table(tab, sep='\s+', skipinitialspace=True)
        data = data.reset_index(drop=True)
        data['chrom'] = data['chrom'].astype(str)
        
        self.compare_tables(severity.head(), data)
