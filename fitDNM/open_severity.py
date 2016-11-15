
import tempfile

import pandas
import pysam

from fitDNM.mutation_rates import flatten_indexed_table

def get_cadd_severity(symbol, chrom, start, end,
        cadd_path='/lustre/scratch113/projects/ddd/users/ps14/CADD/whole_genome_SNVs.tsv.gz'):
    ''' get per nucleotide mutation rates for all SNV alt alleles in a gene
    
    
    See downloadable CADD files here: http://cadd.gs.washington.edu/download
    
    Args:
        symbol: HGNC symbol for gene
        chrom: chromosome
        start: start position
        end: end position
        cadd_path: path to CADD file (bgzip compressed and tabix-indexed)
    
    Returns:
        pandas DataFrame of severity scores at each possible SNV change within
        a genome region.
    '''
    
    tabix = pysam.TabixFile(cadd_path)
    
    with tempfile.TemporaryFile() as handle:
        for x in tabix.fetch(chrom, start, end):
            handle.write(x + '\n')
        
        # ensure all the data is written, and we can read from the start
        handle.flush()
        handle.seek(0)
        
        cadd = pandas.read_table(handle, header=None,
            names=['chrom', 'pos', 'ref', 'alt', 'raw', 'score'])
    
    # scale the cadd scores between 0-1
    cadd['score'] = 1 - 10**-(cadd['score']/10)
    cadd['gene'] = symbol
    
    cadd = cadd.pivot_table(rows=['chrom', 'pos', 'gene', 'ref'],
        cols='alt', values='score' ,fill_value=0.0)
    
    return flatten_indexed_table(cadd)
