
library(argparse)
library(fitDNM)

get_options <- function() {
    parser = ArgumentParser()
    parser$add_argument('--males', type='integer', help='number of males.')
    parser$add_argument('--females', type='integer', help='number of females.')
    parser$add_argument('--de-novos', help='Path to table of de novos.')
    parser$add_argument('--rates', help='Path to table of per-base and allele mutation rates.')
    parser$add_argument('--severity', help='Path to table of per base and allele severity scores.')
    parser$add_argument('--output', help='Path to put output files into.')

    args = parser$parse_args()
    
    return(args)
}

main <- function() {
    args = get_options()
    
    # args$males = 156
    # args$females = 108
    # args$de_novos = '../data-raw/de_novos.txt'
    # args$rates = '../data-raw/rates.txt.gz'
    # args$severity = '../data-raw/severity.txt.gz'
    # args$output = '../output.txt'
    
    de_novos = read.table(args$de_novos, header=TRUE, stringsAsFactors=FALSE)
    mu_rate = read.table(args$rates, header=TRUE, stringsAsFactors=FALSE)
    severity = read.table(args$severity, header=TRUE, stringsAsFactors=FALSE, fill=TRUE)
    
    exclude = get_rows_to_exclude(severity)
    
    mu_rate = tidy_table(mu_rate, exclude)
    severity = tidy_table(severity, exclude)
    
    computed = list()
    for (symbol in sort(unique(de_novos[['gene']]))) {
        values = compute_pvalue(de_novos, args$males, args$females, symbol,
                severity, mu_rate)
        computed[[symbol]] = values
    }
    
    # convert the output to a table and save to disk
    computed = do.call(rbind, computed)
    write.table(computed, file=args$output, sep='\t', quote=FALSE,
        row.names=FALSE)
}

main()
