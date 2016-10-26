
import os
import argparse

import pandas

from fitDNM.prepare_data import get_rows_to_exclude, tidy_table
from fitDNM.gene_enrichment import compute_pvalue
from fitDNM.mutation_rates import get_gene_rates

def get_options():
    ''' parse the command line arguments
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--males', type=int, help='number of males.')
    parser.add_argument('--females', type=int, help='number of females.')
    parser.add_argument('--de-novos', help='Path to table of de novos.')
    parser.add_argument('--severity', help='Path to table of per base and allele severity scores.')
    parser.add_argument('--output', help='Path to put output files into.')
    
    parser.add_argument("--genome-build", dest="genome_build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build " \
        "that the de novo coordinates are based on (GrCh37 or GRCh38")
    parser.add_argument("--cache-folder", dest="cache_dir",
        default=os.path.join(os.path.dirname(__file__), "cache"), help="folder" \
        "to cache Ensembl data into (defaults to clustering code directory)")
    
    return parser.parse_args()

def main():
    args = get_options()
    
    # args.males = 156
    # args.females = 108
    # args.de_novos = '../data-raw/de_novos.txt'
    # args.rates = '../data-raw/rates.txt.gz'
    # args.severity = '../data-raw/severity.txt.gz'
    # args.output = '../output.txt'
    
    de_novos = pandas.read_table(args.de_novos)
    
    severity = pandas.read_table(args.severity, compression='gzip')
    exclude = get_rows_to_exclude(severity)
    severity = tidy_table(severity, exclude)
    
    computed = []
    for symbol in sorted(set(de_novos['gene'])):
        print(symbol)
        
        try:
            mu_rate = get_gene_rates(symbol, de_novos)
        except IndexError:
            continue
        
        values = compute_pvalue(de_novos, args.males, args.females, symbol,
                severity, mu_rate)
        
        if values is not None:
            computed.append(values)
        else:
            print('cannot find mutation rates or severity scores for {}'.format(symbol))
    
    # convert the output to a table and save to disk
    computed = pandas.DataFrame(computed, columns=['symbol', 'cohort_n',
        'nsnv_o', 'n_sites', 'n_de_novos', 'scores', 'p_value', 'p_unweighted'])
    computed.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
