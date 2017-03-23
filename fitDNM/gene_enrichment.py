
import math

from numpy import array
from scipy.stats import poisson

from fitDNM.double_saddle_point import double_saddle_point_approximation

def downweight_severity(data):
    ''' downweight severities by the rate across the alleles at each site
    
    severity = severity - (severity * mu_rate).apply(sum, axis=0)
    '''
    
    data['value'] = data['prob'] * data['score']
    recode = dict([ (pos, sum(x['value'])) for pos, x in data.groupby('pos') ])
    data['score'] = data['score'] - data['pos'].map(recode)
    del data['value']
    
    return data

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
    
    data = data[~data[['pos', 'alt']].duplicated()]
    return data.sort_values(['pos', 'alt'])

def get_expected_rates(de_novos, data, n_male, n_female):
    '''
    '''
    
    # count the samples in the cohort, so we can calculate expected mutations
    # TODO: add in a better chromX adjustment
    nsample = n_male + n_female
    if set(de_novos['chrom']) == set(['X']):
        nsample = n_male / 2 + n_female
    
    # get the expected mutation rates per base per site for the gene
    return nsample * 2 * data['prob']

def compute_pvalue(de_novos, n_male, n_female, symbol, severity, rates):
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
    data = downweight_severity(data)
    data['prob'] = get_expected_rates(de_novos, data, n_male, n_female)
    
    if any(data['score'].isnull()):
        print("{} has NA in severity".format(symbol))
    
    # intersect the de novos with rate and severity data. This restricts
    # de novos to be within the coding sequence of the gene.
    de_novos = de_novos.merge(data, on=['gene', 'chrom', 'pos', 'ref', 'alt'])
    observed_score = sum(de_novos['score'])
    
    p_value = double_saddle_point_approximation(observed_score, data['prob'], data['score'])
    p_unweighted = poisson.sf(len(de_novos) - 1, sum(data['prob']))
    
    return {'symbol': symbol, 'nsnv_o': sum(data['score']),
        'n_sites': len(set(data['pos'])), 'n_de_novos': len(de_novos),
        'scores': round(observed_score, 3), 'p_value': p_value,
        'p_unweighted': p_unweighted}
