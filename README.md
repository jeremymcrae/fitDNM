### fitDNM (python-style)
This repository contains the [fitDNM *de novo* enrichment model](http://dx.doi.org/10.1016/j.ajhg.2015.06.013)
converted to python.

In addition to the python conversion, this code changes some inputs:
 - loads mutation rates on the fly
 - uses CADD severity scores rather than PolyPhen scores.
 - allows for longer sequence context mutation rate models (e.g. 7-mer model)

#### Installation
``` sh
python setup.py install --user
```

#### Example files
 - `examples/de_novos.txt`: an example de novo mutation file for input of fitDNM
 - `examples/output.txt`: an example output file of fitDNM

#### Running
``` sh
python scripts/run_fitDNM.py \
  --males 156 \
  --females 108 \
  --de-novos examples/de_novos.txt \
  --severity CADD_SNV_PATH \
  --output examples/output.txt
```

You can optionally use your own sequence context based rates with
`--rates RATES_PATH`, giving a path to a table with three columns (named 'from',
'to', and 'mu_rate'). 'from' contains the sequence context surrounding the base
to be altered (altered base in the center), 'to' contains the sequence following
the base change and 'mu_rate' contains the numerical mutation rate.

See script documentation with `python scripts/run_fitDNM.py --help`. Alternatively,
the arguments passed to `run_fitDNM.py` are:
 - `--males` number of males in the study. e.g. 156 in example data
 - `--females` number of females in the study. e.g. 108 in example data
 - `--de-novos` file lists de novo mutations information with title line
    (chr pos ref alt gene) Note: we cannot process chromosome Y now. e.g.
    examples/de_novos.txt
 - `--severity` path to table of CADD scores for SNVs throughout the genome. See
    http://cadd.gs.washington.edu/download. The file must be bgzip compressed and
    tabix-indexed for rapid access.
 - `--output` file to output results. e.g. output.txt
