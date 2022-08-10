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
