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
