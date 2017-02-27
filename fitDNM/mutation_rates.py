

import pandas

from denovonear.load_gene import load_gene, \
    get_de_novos_in_transcript
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates

def get_gene_rates(symbol, de_novos, ensembl=None, mut_path=None):
    ''' get per nucleotide mutation rates for all SNV alt alleles in a gene
    
    Args:
        symbol: HGNC symbol for gene
        de_novos: pandas DataFrame containing de novo candidates.
        ensembl: EnsemblRequest object, for extracting gene coordinates and sequence.
        mut_path: path to table of sequence context based mutation rates.
    
    Returns:
        pandas DataFrame of mutation rates at each possible SNV change within
        the coding sequence of a gene.
    '''
    
    mut_dict = load_mutation_rates(mut_path)
    
    if ensembl is None:
        ensembl = EnsemblRequest('cache', 'grch37')
    
    positions = de_novos['pos'][de_novos['gene'] == symbol]
    transcripts = load_gene(ensembl, symbol, positions)
    chrom = transcripts[0].get_chrom()
    
    mu_rate = []
    for transcript in transcripts:
        rates = SiteRates(transcript, mut_dict)
        for cq in ['nonsense', 'missense', 'synonymous', 'splice_lof', 'splice_region']:
            for choice in rates[cq]:
                choice['pos'] = transcript.get_position_on_chrom(choice['pos'], choice['offset'])
                choice['gene'] = symbol
                choice['chrom'] = chrom
                mu_rate.append(choice)
    
    # convert the list of dictionaries to a DataFrame, then reshape to the
    # required form for fitDNM
    mu_rate = pandas.DataFrame(mu_rate)
    mu_rate = mu_rate.pivot_table(index=['gene', 'chrom', 'pos', 'ref'],
        columns='alt', values='prob')
    
    mu_rate = mu_rate.fillna(0)
    
    return flatten_indexed_table(mu_rate)

def flatten_indexed_table(data):
    ''' integrates MultiIndex columns in a DataFrame into the table
    
    Args:
        data: pandas DataFrame, with a MultiIndex following a pivot_table()
            reshaping of a DataFrame.
    
    Returns:
        table where the index columns are now columns within the table
    '''
    
    idx_names = list(data.index.names)
    idx_data = [ list(data.index.get_level_values(x)) for x in idx_names ]
    
    col_names = list(data.columns)
    col_data = [ list(data[x]) for x in col_names ]
    
    # convert to a dictionary, indexed by column names, to use in DataFrame
    flattened = dict(zip(idx_names + col_names, idx_data + col_data))
    
    return pandas.DataFrame(flattened, columns=idx_names + col_names)
