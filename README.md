# SComatic
A tool for detecting somatic variants in single cell data

## Installation and requirements
The code requires python version >=3.7.0 and R version >=3.5.0.

## Step 1: Split alignment file in cell type specific bams

The BAM file containing the sequencing reads for all cell types in a sample are split into cell-type-specific BAM files using precomputed cell type annotations. It requires to have the cell type barcode information in the tag bam file “CB”. In addition, the precomputed cell type annotations should be provided in the next format: 

The command line to run this step is: 

- List of parameters:
```
python SplitBam/SplitBamCellTypes.01102021.py --help
usage: SplitBamCellTypes.01102021.py [-h] --bam BAM --meta META [--id ID]
                                     [--tissue TISSUE] [--max_nM MAX_NM]
                                     [--min_MQ MIN_MQ] [--n_trim N_TRIM]
                                     [--outdir OUTDIR]

Split bam file by cell types

optional arguments:
  -h, --help       show this help message and exit
  --bam BAM        Bam file to be analysed (Sorted by coordinates)
  --meta META      Metadata with cell barcodes per cell type
  --id ID          Sample id
  --tissue TISSUE  Tissue name. Recommended when different tissues from the
                   same individual are analysed
  --max_nM MAX_NM  Maximum number of mismatched permitted to consider the read
                   [Default: 1000]
  --min_MQ MIN_MQ  Minimum mapping quality required to consider the read
                   [Default: 255]
  --n_trim N_TRIM  Number of bases trimmed (setting the base quality to 0) at
                   the beginning and end of the read [Default: 0]
  --outdir OUTDIR  Out directory
```

- Run example:
```python
python SplitBam/SplitBamCellTypes.01102021.py --bam $bam \
        --meta $out1/${sample}.cell_annotation.ready.txt \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --outdir $out2
 ```
 
## Step 2: Collecting base count information

The count of each base in each cell type for every position in the genome is recorded in a base count matrix indexed by cell types and genomic coordinates

The command line to run this step is: 

- List of parameters:
```
python BaseCellCounter/BaseCellCounter.29032021.py --help
usage: BaseCellCounter.29032021.py [-h] --bam BAM --ref REF --chrom CHROM
                                   [--out_folder OUT_FOLDER] [--id ID]
                                   [--nprocs NPROCS] [--bin BIN] [--bed BED]
                                   [--bed_out BED_OUT] [--min_ac MIN_AC]
                                   [--min_af MIN_AF] [--min_dp MIN_DP]
                                   [--min_cc MIN_CC] [--min_bq MIN_BQ]
                                   [--min_mq MIN_MQ] [--tmp_dir TMP_DIR]

Script to obtain a list of base and cell counts in scRNA bam file

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM             Tumor bam file to be analysed
  --ref REF             Reference genome. *fai must be available in the same
                        folder as reference
  --chrom CHROM         Chromosome to be analysed. --chrom all to run
                        in all chromosomes
  --out_folder OUT_FOLDER
                        Out folder
  --id ID               Prefix of out file. If provided, please use next
                        format: *.[cell_type] . Example: sample1.t_cell. If
                        not provided, it is taken from bam file
  --nprocs NPROCS       Number of processes [Default = 1]
  --bin BIN             Bin size for running the analysis [Default 50000]
  --bed BED             Regions to focus the analysis. Three-column bed file
  --bed_out BED_OUT     Regions to ignore in the analysis. Three-column bed
                        file
  --min_ac MIN_AC       Minimum alternative count to consider the genomic
                        site. Default = 0
  --min_af MIN_AF       Minimum alternative allele fraction to consider the
                        genomic site. Default = 0
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

- Run example:
```python
python BaseCellCounter/BaseCellCounter.29032021.py --bam $cell_type_bam \
    --ref $REF \
    --chrom all \
    --bed_out $BLACKLIST \
    --out_folder $out2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 16
```

If running in an HPC cluster, one might want to compute the base count information for all cell types at once. One can do that by running, for example:

```python
for bam in $(ls -d $out/*/cell_type_bams/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Out for bam cell counter
  out2=$(dirname $(dirname $bam))/base_cell_counter
  mkdir -p $out2

  # Log folder
  logs=$out2/logs
  mkdir -p $logs

  # Temp folder
  temp=$out2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  run="python $SCRIPT --bam $bam \
    --ref $REF \
    --chrom all \
    --bed_out $BLACKLIST \
    --out_folder $out2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 16"

  rm -f ${logs}/BaseCellCounter.o.log ${logs}/BaseCellCounter.e.log
  bsub -M 5G -n 16 -J BaseCellCounter -o ${logs}/BaseCellCounter.${cell_type}.o.log -e ${logs}/BaseCellCounter.${cell_type}.e.log "$run"
done
```

## Step 3: Merge all base counts in a multiple cell type matrix
This script takes all computed cell type matrixes and merge them in a unique one (tsv format). All cell type matrixes must be in the same folder, and non-desired tsv files must be located in this folder. 

- List of parameters:
```python 
python MergeCounts/MergeBaseCellCounts.29032021.py --help
usage: MergeBaseCellCounts.29032021.py [-h] --tsv_folder TSV_FOLDER --out_file
                                       OUT_FILE

Script to merge the cell/base counts tsv files per cell type in only one

optional arguments:
  -h, --help            show this help message and exit
  --tsv_folder TSV_FOLDER
                        Folder with cell/base count files (tsv) per cell type.
                        Avoid not desired tsv files in this folder
  --out_file OUT_FILE   Prefix out files
```

- Run example
```python
python MergeCounts/MergeBaseCellCounts.29032021.py --tsv_folder $TSV \
  --out_file $out2/${sample}.BaseCellCounts.AllCellTypes.tsv
```

## Step 4: Detection of somatic variants
The variant calling takes the merged base count matrix to finally provide a final list of somatic variants. It is split in two steps:

### Step 4.1

Firstly, it applies a set of hard filters and the beta-binomial tests to distinguish background error noise from potential somatic mutations. 

- This script takes the next parameters:
```python
python BaseCellCalling/BaseCellCalling.step1.01102021.py --help
usage: BaseCellCalling.step1.01102021.py [-h] --infile INFILE --outfile
                                         OUTFILE --ref REF [--editing EDITING]
                                         [--pon PON] [--min_cov MIN_COV]
                                         [--min_cells MIN_CELLS]
                                         [--min_ac_cells MIN_AC_CELLS]
                                         [--min_ac_reads MIN_AC_READS]
                                         [--max_cell_types MAX_CELL_TYPES]
                                         [--min_cell_types MIN_CELL_TYPES]
                                         [--fisher_cutoff FISHER_CUTOFF]
                                         [--min_distance MIN_DISTANCE]
                                         [--alpha1 ALPHA1] [--beta1 BETA1]
                                         [--alpha2 ALPHA2] [--beta2 BETA2]

Script to perform the scRNA somatic variant calling

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       Input file with all samples merged in a single tsv
  --outfile OUTFILE     Out file prefix
  --ref REF             Reference fasta file (*fai must exist)
  --editing EDITING     Editing file to be used to remove RNA-diting sites
  --pon PON             Panel of normals (PoN) file to be used to remove
                        germline and false positive calls
  --min_cov MIN_COV     Minimum depth of coverage to consider a sample.
                        [Default: 5]
  --min_cells MIN_CELLS
                        Minimum number of cells covering a site to consider a
                        sample. [Default: 5]
  --min_ac_cells MIN_AC_CELLS
                        Minimum number of cells supporting the alternative
                        allele to consider a mutation. [Default: 2]
  --min_ac_reads MIN_AC_READS
                        Minimum number of reads supporting the alternative
                        allele to consider a mutation. [Default: 3]
  --max_cell_types MAX_CELL_TYPES
                        Maximum number of cell types carrying a mutation to
                        make a somatic call. [Default: 1]
  --min_cell_types MIN_CELL_TYPES
                        Minimum number of cell types with enough coverage and
                        cell to consider a site as callable [Default: 2]
  --fisher_cutoff FISHER_CUTOFF
                        P-value cutoff for the Fisher exact test performed to
                        detect strand bias. Expected float value, if applied,
                        we recommend 0.001. By default, this test is switched
                        off with a value of 1 [Default: 1]
  --min_distance MIN_DISTANCE
                        Minimum distance allowed between potential somatic
                        variants [Default: 5]
  --alpha1 ALPHA1       Alpha parameter for beta distribution of read counts.
                        [Default: 0.260288007167716]
  --beta1 BETA1         Beta parameter for beta distribution of read counts.
                        [Default: 173.94711910763732]
  --alpha2 ALPHA2       Alpha parameter for beta distribution of cell counts.
                        [Default: 0.08354121346569514]
  --beta2 BETA2         Beta parameter for beta distribution of cell counts.
                        [Default: 103.47683488327257]
```

- Run example:
```python
python BaseCellCalling/BaseCellCalling.step1.01102021.py \
          --infile $TSV \
          --outfile $out2/${sample} \
          --ref $REF
```

In case that the user wants to estimate new beta-binomial parameters, this extra step should run (LINK).

### Step 4.2 

Secondly, it takes the output of the previous step and applies other filters based on external datasets (RNA editing and Panel of normals), as well as labelling clustered variants. High quality variants will show the label “PASS” in the FILTER column of the output file. All columns and other variables presented in the final file are described in the header of the file (vcf-like file).

- This script has these parameters: 
```python 
python BaseCellCalling/BaseCellCalling.step2.01102021.py --help
usage: BaseCellCalling.step2.01102021.py [-h] --infile INFILE --outfile
                                         OUTFILE [--editing EDITING]
                                         [--pon PON]
                                         [--min_distance MIN_DISTANCE]

Script to perform the scRNA somatic variant calling

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       Input file with all samples merged in a single tsv
  --outfile OUTFILE     Out file prefix
  --editing EDITING     Editing file to be used to remove RNA-diting sites
  --pon PON             Panel of normals (PoN) file to be used to remove
                        germline and false positive calls
  --min_distance MIN_DISTANCE
                        Minimum distance allowed between potential somatic
                        variants [Default: 5]
```

- Run example:
```python
python BaseCellCalling/BaseCellCalling.step2.01102021.py \
          --infile $out2/${sample}.step1.targeted_regions.tsv \
          --outfile $out2/${sample}.targeted_regions  \
          --editing $editing \
          --pon $PON
```

## Other tools of interest

### Estimate new beta-binomial parameters
To estimate the beta-binomial parameters in a new set of samples, take a look at [Beta-binomial estimation](/docs/betabinomialestimation.md) for more information.

### Panel of normals
You can access the panel of normals description here [PON](/docs/pon.md) for more information.

### Get the number of callable sites per cell type

### Get the number of callable sites per unique cell

### Get the genotype for each unique cell in variant sites

### Prepare output file for annovar annotation

