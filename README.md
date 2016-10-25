This directory contain scripts and examples to run fitDNM

Files included in this directory:
 - `data-raw/de_novos.txt`: an example de novo mutation file for input of fitDNM
 - `data-raw/rates.txt.gz`: Example of mutation rate file
 - `data-raw/severity.txt.gz`: Example of severity file

#### Installation
``` sh
python setup.py install --user
```

#### Running
``` sh
python scripts/run_fitDNM.py \
  --males 156 \
  --females 108 \
  --de-novos data-raw/de_novos.txt \
  --rates data-raw/rates.txt.gz \
  --severity data-raw/severity.txt.gz \
  --output data-raw/output.txt
```

See script documentation with `python scripts/run_fitDNM.py --help`. Alternatively,
the arguments passed to `run_fitDNM.py` are:
 - `--males` number of males in the study. e.g. 156 in example data
 - `--females` number of females in the study. e.g. 108 in example
 - `--de-novos` file contains the de novo mutation information with title line
    (chr pos ref alt gene) Note: we cannot process chromosome Y now. e.g.
    data-raw/de_novos.txt
 - `--rates`: mutation rate for each locus. Note, the reference position is labeled
    as mutation rate 0. e.g. data-raw/rates.txt.gz
 - `--severity` The polyphen-2 score for all loci of genes used in analysis.
    (-0.001 refer to Polyphen-2 score is not available, synonymous mutations are
    coded as 3, LOF mutations are coded as 2) e.g. data-raw/severity.txt.gz
 - `--output` file to output results. e.g. output.txt
