# unit testing for the fitDNM functions

import unittest

import pandas

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from fitDNM.prepare_data import get_rows_to_exclude, tidy_table
from tests.compare_dataframes import CompareTables

class TestPrepareDataPy(CompareTables):
    ''' Data cleaning checks
    '''
    
    def test_get_rows_to_exclude(self):
        ''' get_rows_to_exclude is correct when all rows are fine
        '''
        
        tab = StringIO('''gene  chr  pos   A  T  G  C
                         GENE1   1    1    1  0  0  0
                         GENE1   1    2    1  0  0  0
                         GENE1   1    2    1  0  0  0''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        
        self.assertEqual(list(get_rows_to_exclude(data)), [False, False, False])
        
    def test_get_rows_to_exclude_some_bad(self):
        ''' get_rows_to_exclude is correct when all some rows have bad values
        '''
        tab = StringIO('''gene  chr  pos   A  T  G   C
                         GENE1   1    1   -1  0  0   0
                         GENE1   1    2    1  0  0   0
                         GENE1   1    2    1  0  NA  0''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        
        self.assertEqual(list(get_rows_to_exclude(data)), [True, False, True])
    
    def test_tidy_table(self):
        ''' check that tidy_table is correct
        '''
        
        tab = StringIO('''Gene  Chr  Pos   Ref  extra
                           GENE1   1    1    A     0
                           GENE1   1    2    T     0
                           GENE1   1    2    G     0''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        data = data.reset_index(drop=True)
        exclude = pandas.Series([False, False, False])
        
        expected = data.rename(columns={'Gene': 'gene', 'Chr': 'chr',
            'Pos': 'pos', 'Ref': 'ref'})
        
        self.compare_tables(tidy_table(data, exclude), expected)
    
    def test_tidy_table_exclude_rows(self):
        ''' check that tidy_table is correct when we have rows to exclude
        '''
        
        tab = StringIO('''Gene  Chr  Pos   Ref  extra
                           GENE1   1    1    A     0
                           GENE1   1    2    T     0
                           GENE1   1    2    G     0''')
        
        data = pandas.read_table(tab, delim_whitespace=True, skipinitialspace=True)
        data = data.reset_index(drop=True)
        exclude = pandas.Series([False, True, False])
        
        expected = data.rename(columns={'Gene': 'gene', 'Chr': 'chr',
            'Pos': 'pos', 'Ref': 'ref'})
        expected = expected.iloc[[0, 2]]
        
        self.compare_tables(tidy_table(data, exclude), expected)
