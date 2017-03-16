
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
    
    def parse(line):
        chrom, pos, ref, alt, raw, scaled = line.split('\t')
        return {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt,
            'raw': float(raw), 'score': float(scaled)}
        
    cadd = pandas.DataFrame([ parse(x) for x in tabix.fetch(chrom, start, end) ])
    
    # scale the cadd scores between 0-1
    cadd['score'] = 1 - 10**-(cadd['score']/10)
    cadd['gene'] = symbol
    cadd['chrom'] = cadd['chrom'].astype(str)
    
    cadd = cadd.pivot_table(index=['chrom', 'pos', 'gene', 'ref'],
        columns='alt', values='score' ,fill_value=0.0)
    
    return flatten_indexed_table(cadd)
