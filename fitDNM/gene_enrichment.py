
import math

from numpy import array
from scipy.stats import poisson

from fitDNM.saddlepoint import saddlepoint

def downweight_severity(data):
    ''' downweight severities by the rate across the alleles at each site
    
    severity = severity - (severity * mu_rate).apply(sum, axis=0)
    '''
    
    data = data.copy()
    data['value'] = data['prob'] * data['score']
    recode = dict([ (pos, sum(x['value'])) for pos, x in data.groupby('pos') ])
    data['score'] = data['score'] - data['pos'].map(recode)
    del data['value']
    
    return data['score']

def harmonise_data(rates, severity, symbol):
    ''' merge mutation rates and severity score datasets
    
    Args:
        rates: pandas DataFrame that includes a 'prob' column
        severity: pandas DataFrame that includes a 'score' column
        symbol: HGNC symbol for gene
    
    Returns:
        tuple of (mutation probabilities, severity score, position, )
    '''
    
    if len(rates) == 0 or len(severity) == 0:
        raise ValueError('empty rates or empty severity dataset!')
    
    data = rates.merge(severity, on=['gene', 'chrom', 'pos', 'ref', 'alt'])
    
    if len(data) == 0:
        raise ValueError('no shared sites between rates and severity datasets!')
    
    data = data[data['gene'] == symbol]
    
    assert len(set(data['chrom'])) == 1
    
    data = data[~data[['pos', 'alt']].duplicated()]
    return data.sort_values(['pos', 'alt'])

def get_expected_rates(data, male, female):
    ''' calculate expected number of mutations per site.
    
    Args:
        data: pandas DataFrame, which includes columns for 'chrom' (for
            chromosome) and 'prob' (for mutation probability)
        male: number of male probands
        female: number of female probands
    
    Returns:
        Mutation probabilies scaled so that each is the chance of observing a
        mutation at that site, given the number of male and female probands in
        the cohort.
    '''
    
    # count the samples in the cohort, so we can calculate expected mutations
    # TODO: add in a better chromX adjustment
    count = male + female
    if set(data['chrom']) == set(['X']):
        count = male / 2 + female
    
    # get the expected mutation rates per base per site for the gene
    return count * 2 * data['prob']

def enrichment(de_novos, n_male, n_female, symbol, severity, rates):
    ''' compute de novo enrichment for a gene
    
    Args:
        de_novos: table of de novos
        n_male: number of males
        n_female: number of females
        symbol: HGNC symbol for a gene
        severity: table of per base and per allele severity scores, for every gene
        rates: table of mutation rates
    
    Returns:
        dictionary of de novo enrichment results for a gene
    '''
    
    data = harmonise_data(rates, severity, symbol)
    
    data['score'] = downweight_severity(data)
    data['prob'] = get_expected_rates(data, n_male, n_female)
    
    if any(data['score'].isnull()):
        print("{} has NA in severity".format(symbol))
    
    # intersect the de novos with rate and severity data. This restricts
    # de novos to be within the coding sequence of the gene.
    data_prefixed = data.chrom.loc[0].startswith('chr')
    dnm_prefixed = de_novos.chrom.loc[0].startswith('chr')
    if data_prefixed and not dnm_prefixed:
        de_novos['chrom'] = 'chr' + de_novos['chrom']
    de_novos = de_novos.merge(data, on=['gene', 'chrom', 'pos', 'ref', 'alt'])
    observed_score = sum(de_novos['score'])
    
    p_value = saddlepoint(observed_score, data['prob'], data['score'])
    p_unweighted = poisson.sf(len(de_novos) - 1, sum(data['prob']))
    
    return {'symbol': symbol, 'gene_scores': sum(data['score']),
        'sites': len(set(data['pos'])), 'de_novos': len(de_novos),
        'de_novos_score': round(observed_score, 3), 'p_value': p_value,
        'p_unweighted': p_unweighted}
