
import pandas

def get_gene_coords(symbol, gencode):
    ''' get the genomic coordinates for a HGNC symbol from gencode
    
    Args:
        symbol: HGNC symbol as string
        gencode: pandas DataFrame of gencode coordinates (with 'gene', 'chrom',
            'start' and 'end' columns)
    
    Returns:
        tuple of chrom, star and end coordinates
    '''
    
    # TODO: what happens when symbol isn't in gencode, or present more than once?
    row = gencode[gencode['gene'] == symbol].squeeze()
    
    return row['chrom'], row['start'], row['end']

def open_gencode(path):
    ''' get gene ranges from table of gencode coordinates
    
    Args:
        path: path to gencode table
    
    Returns:
        pandas DataFrame of gencode coordinates, restricted to protein coding genes
    '''
    
    path = '/lustre/scratch113/projects/ddd/users/jm33/all_gencode_genes_v19.txt.gz'
    
    gencode = pandas.read_table(path, compression='gzip')
    
    gencode = gencode[gencode['gene_type'].isin(['protein_coding', 'polymorphic_pseudogene'])]
    
    # strip the 'chr' from the chromosome string to match other inputs
    gencode['chrom'] = gencode['chr'].str.replace('chr', '')
    gencode['end'] = gencode['stop']
    
    return gencode[['gene', 'chrom', 'start', 'end']].copy()
