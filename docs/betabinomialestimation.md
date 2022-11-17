## Estimating parameters for the beta-binomial distribution

To estimate the Beta-binomial distribution parameters in new data, the user needs to provide (--in_tsv) a tsv (or txt) file listing the full path to all desired sample-cell-type (ideally cancer-free samples) obtained by BaseCellCounter.py. The file should look like next:

```
/path/to/sample1.epithelial.tsv
/path/to/sample1.T_cells.step1.tsv
/path/to/sample1.myeloid.step1.tsv
/path/to/sample2.epithelial.tsv
/path/to/sample2.T_cells.step1.tsv
/path/to/sample2.myeloid.step1.tsv
/path/to/sample3.epithelial.tsv
/path/to/sample3.T_cells.step1.tsv
/path/to/sample3.myeloid.step1.tsv
```

The python script takes these parameters:


```
python scripts/BetaBinEstimation/BetaBinEstimation.py --help
usage: BetaBinEstimation.py [-h] --in_tsv IN_TSV --outfile OUTFILE
                            [--n_sites N_SITES] [--seed SEED]

Script to estimate the Beta-binomial distribution parameters (alpha and beta)
to be used afterwards in the BaseCellCalling.step1.py

optional arguments:
  -h, --help         show this help message and exit
  --in_tsv IN_TSV    File listing the tsv files to be used for the beta-
                     binomial fitting (obtained with BaseCellCounter.py
                     script)
  --outfile OUTFILE  Report with the estimated Beta-binomial parameters
  --n_sites N_SITES  Approximate number of sites to be used for fitting the
                     Beta-binomial distribution [Default: 500000]
  --seed SEED        Random seed for computation [Default: 1992]
```
