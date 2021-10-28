
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

def weight_site(site, weights, lof_weight):
    ''' assign a weight to the the severity score
    '''
    
    cq = site['consequence']
    
    # check whether the site is under regional constraint
    constraint = 'unconstrained'
    if site['constrained']:
        constraint = 'constrained'
    
    if cq == 'nonsense' or cq == 'splice_lof':
        return lof_weight
    else:
        for (low, high) in weights[constraint]:
            if low <= site['score'] < high:
                return weights[constraint][(low, high)]
    
    return None

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
    
    weights = {
        'unconstrained': {
            (0, 2): 0.8937899916955947, (2, 4): 1.2018380935433957,
            (4, 6): 1.9210019744290012, (6, 8): 2.095559182524437,
            (8, 10): 1.8958839116730013, (10, 12): 1.2379790547896217,
            (12, 14): 1.1536797961075576, (14, 16): 0.9571549689252671,
            (16, 18): 1.0637563686850184, (18, 20): 1.5602804818006433,
            (20, 22): 2.1320822607373815, (22, 24): 2.8802645134215537,
            (24, 26): 4.490358805819542, (26, 28): 5.0898990323887325,
            (28, 30): 6.248469361168584, (30, 32): 6.836303061913207,
            (32, 34): 6.367073811801687, (34, 36): 5.687929107769973,
            (36, 60): 5.687929107769973},
        'constrained': {
            (0, 2): 0.000000000000000, (2, 4): 8.182502325314887,
            (4, 6): 5.223168595924202, (6, 8): 0.000000000000000,
            (8, 10): 1.9755254351570928, (10, 12): 2.757582326295284,
            (12, 14): 1.127500058867108, (14, 16): 5.037105250368409,
            (16, 18): 4.713287802915507, (18, 20): 4.480566806148068,
            (20, 22): 7.413052297435974, (22, 24): 9.943230401766982,
            (24, 26): 16.459396897564233, (26, 28): 17.523439447046083,
            (28, 30): 16.7237728797649, (30, 32): 19.445687203439924,
            (32, 34): 24.385071725281048, (34, 36): 35.76543858738138,
            (36, 60): 35.76543858738138}
        }
    
    lof_weight = 30.498635973896338
    
    data = harmonise_data(rates, severity, symbol)
    
    data['score'] = [ weight_site(x, weights, lof_weight) for i, x in data.iterrows() ]
    data['score'] = downweight_severity(data)
    data['prob'] = get_expected_rates(data, n_male, n_female)
    
    if any(data['score'].isnull()):
        print("{} has NA in severity".format(symbol))
    
    # intersect the de novos with rate and severity data. This restricts
    # de novos to be within the coding sequence of the gene.
    de_novos = de_novos.merge(data, on=['gene', 'chrom', 'pos', 'ref', 'alt'])
    observed_score = sum(de_novos['score'])
    
    p_value = saddlepoint(observed_score, data['prob'], data['score'])
    p_unweighted = poisson.sf(len(de_novos) - 1, sum(data['prob']))
    
    return {'symbol': symbol, 'gene_scores': sum(data['score']),
        'sites': len(set(data['pos'])), 'de_novos': len(de_novos),
        'de_novos_score': round(observed_score, 3), 'p_value': p_value,
        'p_unweighted': p_unweighted}
