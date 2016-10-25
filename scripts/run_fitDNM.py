
import argparse

import pandas

from fitDNM.prepare_data import get_rows_to_exclude, tidy_table
from fitDNM.gene_enrichment import compute_pvalue

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument('--males', type=int, help='number of males.')
    parser.add_argument('--females', type=int, help='number of females.')
    parser.add_argument('--de-novos', help='Path to table of de novos.')
    parser.add_argument('--rates', help='Path to table of per-base and allele mutation rates.')
    parser.add_argument('--severity', help='Path to table of per base and allele severity scores.')
    parser.add_argument('--output', help='Path to put output files into.')
    
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
    mu_rate = pandas.read_table(args.rates, compression='gzip')
    severity = pandas.read_table(args.severity, compression='gzip')
    
    exclude = get_rows_to_exclude(severity)
    
    mu_rate = tidy_table(mu_rate, exclude)
    severity = tidy_table(severity, exclude)
    
    computed = []
    for symbol in sorted(set(de_novos['gene'])):
        values = compute_pvalue(de_novos, args.males, args.females, symbol,
                severity, mu_rate)
    
    # convert the output to a table and save to disk
    computed = pandas.DataFrame(computed)
    computed.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
