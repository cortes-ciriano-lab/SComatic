# Other SComatic functionalities

## Computing the number of callable sites per cell type

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

**Example:** check [here](/docs/SComaticExample.md) to see how to run this script with an example sample.  

## Computing the number of callable sites per single cell

```
python scripts/SitesPerCell/SitesPerCell.py --help
usage: SitesPerCell.py [-h] --bam BAM --ref REF [--infile INFILE]
                       [--min_ct1 MIN_CT1] [--min_ct2 MIN_CT2]
                       [--out_folder OUT_FOLDER] [--id ID] [--nprocs NPROCS]
                       [--bin BIN] [--min_dp MIN_DP] [--min_cc MIN_CC]
                       [--min_bq MIN_BQ] [--min_mq MIN_MQ] [--tmp_dir TMP_DIR]

Script to calculate the number of callable sites per unique cell

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Tumor bam file to be analysed
  --ref REF             Reference genome. *fai must be available in the same
                        folder as reference
  --infile INFILE       Base calling file (obtained by
                        BaseCellCalling.step1.py)
  --min_ct1 MIN_CT1     Minimum number of cell types with enough reads to
                        consider a genomic site. Default = 2
  --min_ct2 MIN_CT2     Minimum number of cell types with enough unique cell
                        counts to consider a genomic site. Default = 2
  --out_folder OUT_FOLDER
                        Out folder
  --id ID               Prefix of out file. If provided, please use next
                        format: *.[cell_type] . Example: sample1.t_cell. If
                        not provided, it is taken from bam file
  --nprocs NPROCS       Number of processes [Default = 1]
  --bin BIN             Bin size for running the analysis [Default 50000]
  --min_dp MIN_DP       Minimum coverage to consider the genomic site. Default
                        = 5
  --min_cc MIN_CC       Minimum number of cells required to consider a genomic
                        site. Default = 5
  --min_bq MIN_BQ       Minimum base quality permited for the base counts.
                        Default = 20
  --min_mq MIN_MQ       Minimum mapping quality required to analyse read.
                        Default = 255
  --tmp_dir TMP_DIR     Temporary folder for tmp files
```

**Example:** check [here](/docs/SComaticExample.md) to see how to run this script with an example sample.  

## Computing the genotype for each cell at the variant sites

```
python scripts/SingleCellGenotype/SingleCellGenotype.py --help
usage: SingleCellGenotype.py [-h] --bam BAM --infile INFILE --ref REF --meta
                             META [--out_file OUT_FILE] [--alt_flag {Alt,All}]
                             [--nprocs NPROCS] [--bin BIN] [--min_bq MIN_BQ]
                             [--min_mq MIN_MQ] [--tissue TISSUE]
                             [--tmp_dir TMP_DIR]

Script to get the alleles observed in each unique cell for the variant sites

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Tumor bam file to be analysed
  --infile INFILE       Base calling file (obtained by
                        BaseCellCalling.step2.py), ideally only the PASS
                        variants
  --ref REF             Reference genome. *fai must be available in the same
                        folder as reference
  --meta META           Metadata with cell barcodes per cell type
  --out_file OUT_FILE   Out file
  --alt_flag {Alt,All}  Flag to search for cells carrying the expected alt
                        variant (Alt) or all cells independent of the alt
                        allele observed (All)
  --nprocs NPROCS       Number of processes [Default = 1]
  --bin BIN             Bin size for running the analysis [Default 50000]
  --min_bq MIN_BQ       Minimum base quality permited for the base counts.
                        Default = 30
  --min_mq MIN_MQ       Minimum mapping quality required to analyse read.
                        Default = 255
  --tissue TISSUE       Tissue of the sample
  --tmp_dir TMP_DIR     Temporary folder for tmp files
```

**Example:** check [here](/docs/SComaticExample.md) to see how to run this script with an example sample.  

## Computing the trinucleotide context background

```
python scripts/TrinucleotideBackground/TrinucleotideContextBackground.py --help
usage: TrinucleotideContextBackground.py [-h] --in_tsv IN_TSV --out_file
                                         OUT_FILE

Script to obtain trincucleotide context background

optional arguments:
  -h, --help           show this help message and exit
  --in_tsv IN_TSV      File listing the tsv files to be used for the
                       trinucleotide context background computation (files
                       obtained in BaseCellCalling.step1.py)
  --out_file OUT_FILE  Output file

```

