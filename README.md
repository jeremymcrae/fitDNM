This directory contain scripts and examples to run fitDNM

Files included in this directory:
 - `examples/de_novos.txt`: an example de novo mutation file for input of fitDNM
 - `examples/output.txt`: an example output file of fitDNM

#### Installation
``` sh
python setup.py install --user
```

#### Running
``` sh
python scripts/run_fitDNM.py \
  --males 156 \
  --females 108 \
  --de-novos examples/de_novos.txt \
  --severity CADD_SNV_PATH \
  --output examples/output.txt
```

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
