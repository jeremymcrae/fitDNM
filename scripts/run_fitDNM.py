
import os
import argparse

import pandas

from fitDNM.gene_enrichment import compute_pvalue
from fitDNM.mutation_rates import get_gene_rates
from fitDNM.open_severity import get_cadd_severity

def get_options():
    ''' parse the command line arguments
    '''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--males', type=int, help='number of males.')
    parser.add_argument('--females', type=int, help='number of females.')
    parser.add_argument('--de-novos', help='Path to table of de novos.')
    parser.add_argument('--severity', help='Path to table of per base and allele CADD severity scores.')
    parser.add_argument('--rates', help='Path to table of sequence context '
        'based rates. Defaults to Kaitlin Samocha\'s trinucleotide-based rates.')
    parser.add_argument('--constraint', help='Path to table of regional constraint.')
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
    
    de_novos = pandas.read_table(args.de_novos)
    
    computed = []
    for symbol in sorted(set(de_novos['gene'])):
        print(symbol)
        
        try:
            mu_rate = get_gene_rates(symbol, de_novos, args.constraint, mut_path=args.rates)
        except IndexError:
            continue
        
        chrom, start, end = str(mu_rate['chrom'][0]), min(mu_rate['pos']), max(mu_rate['pos'])
        
        severity = get_cadd_severity(symbol, chrom, start, end, args.severity)
        values = compute_pvalue(de_novos, args.males, args.females, symbol,
                severity, mu_rate)
        computed.append(values)
    
    # convert the output to a table and save to disk
    computed = pandas.DataFrame(computed, columns=['symbol', 'gene_scores',
        'sites', 'de_novos', 'de_novos_score', 'p_value', 'p_unweighted'])
    computed.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()
