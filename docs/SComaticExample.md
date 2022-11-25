# Example of how to run SComatic
This document shows an example of how to run SComatic in a scRNA-seq sample (10x). The provided bam file expands a small fraction of human chromosome 10 (Hg38) and harbours 5 somatic mutations (5 SNVs). The scripts used here are all available at the SComatic/scripts folder, and the example files are located in SComatic/example_data. 

## Prepare files and environment

Activate conda environment if needed
```
conda activate SComatic
```

Create an output folder and go to the main SComatic folder
```
SCOMATIC=SComatic
output_dir=$SCOMATIC/example_data/results
mkdir -p $output_dir
```

If reference genome is not unpacked and indexed (.fai), you have to do it before running SComatic
```
gunzip $SCOMATIC/example_data/chr10.fa.gz
samtools faidx $SCOMATIC/example_data/chr10.fa
```

## Step 1: Splitting alignment file in cell type specific bams
```
sample=Example
output_dir1=$output_dir/Step1_BamCellTypes
mkdir -p $output_dir1

python $SCOMATIC/scripts/SplitBam/SplitBamCellTypes.py --bam $SCOMATIC/example_data/Example.scrnaseq.bam \
        --meta $SCOMATIC/example_data/Example.cell_barcode_annotations.tsv \
        --id ${sample} \
        --n_trim 5 \
        --max_nM 5 \
        --outdir $output_dir1
```

## Step 2: Collecting base count information
```
REF=$SCOMATIC/example_data/chr10.fa

output_dir2=$output_dir/Step2_BaseCellCounts
mkdir -p $output_dir2

for bam in $(ls -d $output_dir1/*bam);do
  
  # Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')

  # Temp folder
  temp=$output_dir2/temp_${cell_type}
  mkdir -p $temp

  # Command line to submit to cluster
  python $SCOMATIC/scripts/BaseCellCounter/BaseCellCounter.py --bam $bam \
    --ref $REF \
    --chrom all \
    --out_folder $output_dir2 \
    --min_bq 30 \
    --tmp_dir $temp \
    --nprocs 1

  rm -rf $temp
done
```

In our experience, when running SComatic in an HPC cluster it is most efficient to compute the base count information for all cell types at once, speacially when submitting each python line independently and at the same time to the cluster, and using multiple processors (p.e. --nprocs 16) 
