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

If running in an HPC cluster, one might want to compute the base count information at once. One can do that by running, for example:

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


