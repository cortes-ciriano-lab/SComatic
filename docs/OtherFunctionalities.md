Getting the number of callable sites per cell type

```
python scripts/GetCallableSites/GetAllCallableSites.py --help
usage: GetAllCallableSites.py [-h] --infile INFILE --outfile OUTFILE
                              [--max_cov MAX_COV]
                              [--min_cell_types MIN_CELL_TYPES]

Script to calculate the number of callable sites per cell type

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       Tsv file obtained by BaseCellCalling.step1.py to be
                        analysed
  --outfile OUTFILE     Out file prefix
  --max_cov MAX_COV     Maximum coverage to record in the callable sites
                        table. Greater values will be collapsed to the
                        provided one. [Default: 150]
  --min_cell_types MIN_CELL_TYPES
                        Minimum number of cell types with enough
                        coverage/cells at a site to be considered as a
                        callable [Default: 2]
```

Getting the number of callable sites per unique cell

Getting the genotype for each unique cell for the variant sites

Preparing the output file for annovar annotation
