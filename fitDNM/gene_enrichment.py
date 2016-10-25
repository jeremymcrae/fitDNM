
from numpy import array
from scipy.stats import poisson

from fitDNM.double_saddle_point import double_saddle_point_approximation

def exclude_and_get_vector(data, exclude):
    ''' exclude certain values from a table ,and convert to a vector
    
    Args:
        data: table of values
        exclude: matching table of booleans, where True indicates values to exclude
    
    Returns:
        list of values, with missing values removed
    '''
    data[exclude] = None
    data = [ x for sublist in data.values.tolist() for x in sublist ]
    
    return array([ x for x in data if x is not None ])

def compute_pvalue(de_novos, n_male, n_female, symbol, severity, mu_rate):
    ''' compute de novo enrichment for a gene
    
    Args:
        de_novos: table of de novos
        n_male: number of males
        n_female: number of females
        symbol: HGNC symbol for a gene
        severity: table of per base and per allele severity scores, for every gene
        mu_rate: table of mutation rates
    
    Returns:
        dictionary of de novo enrichment results for a gene
    '''
    
    # select the mutations for the given gene
    de_novos = de_novos[de_novos['gene'] == symbol]
    
    # count the samples in the cohort, so we can calculate expected mutations
    nsample = n_male + n_female
    if set(de_novos['chr']) == set(['X']):
        nsample = n_male / 2 + n_female
    
    bases = ['A', 'C', 'G', 'T']
    
    severity = severity[severity['gene'] == symbol]
    positions = severity['pos']
    ref_allele = severity['ref']
    severity = severity[bases]
    nsnv_o = severity.sum().sum()
    mu_rate = mu_rate[bases][mu_rate['gene'] == symbol]
    
    # recode LOF (coded as 2) and synonymous (coded as 3) mutations
    severity[severity == 2] = 1
    severity[severity == 3] = 0
    
    if len(severity) == 0:
        return None
    
    # downweight the severities by the rate across the alleles at each site
    severity = severity - (severity * mu_rate).apply(sum, axis=0)
    
    # get the expected mutation rates per base per site for the gene
    mu_rate = nsample * 2 * mu_rate
    
    if any(severity.isnull().apply(any, axis=1)):
        print("{} has NA in severity".format(symbol))
    
    # check if the position is in the CCDS (not included promoter region)
    de_novos = de_novos[de_novos['pos'].isin(positions)]
    
    nsnv = len(severity)
    ndenovo = len(de_novos)
    scores = 0
    
    # if there aren't any sites with mutation rates, or we don't have any de
    # novos, give the gene a p-value of 1
    p_value = 1
    p_unweighted = 1
    if ndenovo > 0:
        for i, row in de_novos.iterrows():
            alt = row['alt']
            if row['ref'] != ref_allele[positions == row['pos']].squeeze():
                raise ValueError('{} {} {} does not match ref'.format(
                    row['chr'], row['pos'],  row['ref']))
            
            scores += severity[alt][positions == row['pos']].squeeze()
        
        # exclude data from non-mutated, i.e. unchanged sites
        exclude = mu_rate == 0
        mu_rate = exclude_and_get_vector(mu_rate, exclude)
        severity = exclude_and_get_vector(severity, exclude)
        
        # compute p-values, using saddlepoint and standard poisson approaches
        p_value = double_saddle_point_approximation(scores, mu_rate, severity)
        p_unweighted = 1 - poisson.pmf(len(de_novos) - 1, sum(mu_rate))
    
    values = {'symbol': symbol, 'cohort_n': nsample, 'nsnv_o': nsnv_o,
        'n_sites': nsnv, 'n_de_novos': ndenovo, 'scores': round(scores, 3),
        'p_value': p_value, 'p_unweighted': p_unweighted}
    
    return values
