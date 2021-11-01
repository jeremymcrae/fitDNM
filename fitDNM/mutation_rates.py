
import asyncio
import logging

import pandas

from denovonear.load_gene import load_gene, \
    construct_gene_object, minimise_transcripts
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.site_specific_rates import SiteRates
from denovonear.rate_limiter import RateLimiter

def get_gene_rates(symbol, de_novos, gencode=None, constraint=None, mut_path=None):
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
        symbol, de_novos, gencode, constraint, mut_path))
    return result

async def _async_get_gene_rates(symbol, de_novos, gencode=None, constraint=None, mut_path=None):
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
    
    positions = de_novos['pos'][de_novos['gene'] == symbol]
    if gencode:
        if symbol not in gencode:
            logging.info(f'cannot find {symbol} in gencode genes')
            raise IndexError
        gene = gencode[symbol]
    else:
        async with RateLimiter(per_second=15) as ensembl:
            gene = await load_gene(ensembl, symbol)
    
    if len(gene.transcripts) == 0:
        logging.info(f'cannot find transcripts for {symbol}')
        raise IndexError
    
    minimized = minimise_transcripts(gene.transcripts, positions)
    transcripts = [x for x in gene.transcripts if x.get_name() in minimized]
    
    if len(transcripts) == 0:
        logging.error(f'DNMs for {symbol} not in any suitable transcript')
        raise IndexError
    
    merged = None
    mu_rate = []
    for transcript in transcripts:
        rates = SiteRates(transcript, mut_dict, merged)
        for cq in ['nonsense', 'missense', 'synonymous', 'splice_lof', 'splice_region']:
            for choice in rates[cq]:
                choice['pos'] = transcript.get_position_on_chrom(choice['pos'], choice['offset'])
                choice['consequence'] = cq
                choice['gene'] = symbol
                choice['chrom'] = gene.chrom
                
                mu_rate.append(choice)
        if merged is None:
            merged = transcript
        merged += transcript
    
    data = pandas.DataFrame(mu_rate)
    data = data.sort_values(['chrom', 'pos', 'alt'])
    data = data.reset_index(drop=True)
    data['chrom'] = data['chrom'].astype(str)
    del data['offset']
    
    return data
