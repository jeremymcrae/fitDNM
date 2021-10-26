
import asyncio

import pandas

from fitDNM.regional_constraint import load_regional_constraint, get_constrained_positions

from denovonear.load_gene import load_gene, \
    construct_gene_object, minimise_transcripts
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.rate_limiter import RateLimiter

def get_gene_rates(symbol, de_novos, constraint=None, mut_path=None):
    ''' get per nucleotide mutation rates for all SNV alt alleles in a gene

    Args:
        symbol: HGNC symbol for gene
        de_novos: pandas DataFrame containing de novo candidates.
        mut_path: path to table of sequence context based mutation rates.

    Returns:
        pandas DataFrame of mutation rates at each possible SNV change within
        the coding sequence of a gene.
    '''
    loop = asyncio.get_event_loop()
    result = loop.run_until_complete(_async_get_gene_rates(
        symbol, de_novos, constraint, mut_path))
    return result

async def _async_get_gene_rates(symbol, de_novos, constraint=None, mut_path=None):
    ''' get per nucleotide mutation rates for all SNV alt alleles in a gene
    
    Args:
        symbol: HGNC symbol for gene
        de_novos: pandas DataFrame containing de novo candidates.
        mut_path: path to table of sequence context based mutation rates.
    
    Returns:
        pandas DataFrame of mutation rates at each possible SNV change within
        the coding sequence of a gene.
    '''
    
    mut_dict = load_mutation_rates(mut_path)
    
    if constraint is not None:
        constraint = load_regional_constraint(constraint)
    else:
        constraint = {'gene': set([])}
    
    async with RateLimiter(per_second=15) as ensembl:
        positions = de_novos['pos'][de_novos['gene'] == symbol]
        gene = await load_gene(ensembl, symbol)
        chrom = gene.chrom
        minimized = minimise_transcripts(gene.transcripts, positions)
        transcripts = [x for x in gene.transcripts if x.get_name() in minimized]
        
        if symbol not in set(constraint['gene']):
            sites = set([])
        else:
            regional = constraint[constraint['gene'] == symbol]
            tx_id = list(regional['transcript'])[0]
            tx = await construct_gene_object(ensembl, tx_id.split('.')[0])
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
