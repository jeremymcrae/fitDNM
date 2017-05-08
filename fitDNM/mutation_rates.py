

import pandas

from fitDNM.regional_constraint import load_regional_constraint, get_constrained_positions

from denovonear.load_gene import load_gene, \
    get_de_novos_in_transcript, construct_gene_object
from denovonear.ensembl_requester import EnsemblRequest
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates

def get_gene_rates(symbol, de_novos, constraint_path, ensembl=None, mut_path=None):
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
    
    constraint = load_regional_constraint(constraint_path)
    positions = de_novos['pos'][de_novos['gene'] == symbol]
    transcripts = load_gene(ensembl, symbol, positions)
    chrom = transcripts[0].get_chrom()
    
    if symbol not in set(constraint['gene']):
        sites = set([])
    else:
        regional = constraint[constraint['gene'] == symbol]
        tx_id = list(regional['transcript'])[0]
        tx = construct_gene_object(ensembl, tx_id.split('.')[0])
        sites = get_constrained_positions(tx, regional, threshold=1e-3, ratio_threshold=0.4)
    
    mu_rate = []
    for transcript in transcripts:
        rates = SiteRates(transcript, mut_dict)
        for cq in ['nonsense', 'missense', 'synonymous', 'splice_lof', 'splice_region']:
            for choice in rates[cq]:
                choice['pos'] = transcript.get_position_on_chrom(choice['pos'], choice['offset'])
                choice['consequence'] = cq
                choice['gene'] = symbol
                choice['chrom'] = chrom
                
                choice['constrained'] = False
                if choice['pos'] in sites:
                    choice['constrained'] = True
                
                mu_rate.append(choice)
    
    data = pandas.DataFrame(mu_rate)
    data = data.sort_values(['chrom', 'pos', 'alt'])
    data = data.reset_index(drop=True)
    data['chrom'] = data['chrom'].astype(str)
    del data['offset']
    
    return data
