docker build -t scomatic .

# Common variables

sample=Example
REF=/opt/SComatic/example_data/chr10.fa
editing=/opt/SComatic/RNAediting/AllEditingSites.hg38.txt
PON=/opt/SComatic/PoNs/PoN.scRNAseq.hg38.tsv
BOI=/opt/SComatic/bed_files_of_interest/UCSC.k100_umap.without.repeatmasker.bed
META=/opt/SComatic/example_data/Example.cell_barcode_annotations.tsv

# Create output directories

output_dir=/home/ubuntu/SComatic/test2/results
mkdir -p $output_dir

# Step 1: Splitting alignment file in cell type specific bams

step1=Step1_BamCellTypes
mount_outdir=/data
output_dir1=$output_dir/$step1
mkdir -p $output_dir1

docker run -v $output_dir:/data -it scomatic SplitBamCellTypes.py \
	--bam /opt/SComatic/example_data/Example.scrnaseq.bam \
	--meta /opt/SComatic/example_data/Example.cell_barcode_annotations.tsv \
	--id ${sample} \
	--n_trim 5 \
	--max_nM 5 \
	--max_NH 1 \
	--outdir $mount_outdir/$step1

# Step 2: Collecting base count information

step2=Step2_BaseCellCounts
mkdir -p $output_dir/$step2

for bam in $(ls -d $output_dir/$step1/*bam);do

	# Cell type
  cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
	echo $cell_type

  # Temp folder
  temp=$output_dir/temp_${cell_type}
  mkdir -p $temp

	fname=$(basename $bam)
	mount_bam="$mount_outdir/$step1/"$fname

	echo $mount_bam
	
  # Command line to submit to cluster
	docker run -v $output_dir:$mount_outdir -it scomatic BaseCellCounter.py --bam $mount_bam \
    --ref $REF \
    --chrom all \
    --out_folder $mount_outdir/$step2  \
    --min_bq 30 \
    --tmp_dir $mount_outdir/temp_${cell_type} \
    --nprocs 1

  rm -rf $temp
done

# Step 3: Merging base count matrices

step3=Step3_BaseCellCountsMerged
output_dir3=$output_dir/$step3
mkdir -p $output_dir3

docker run -v $output_dir:$mount_outdir -it scomatic MergeBaseCellCounts.py \
	--tsv_folder $mount_outdir/$step2 \
  --outfile $mount_outdir/$step3/${sample}.BaseCellCounts.AllCellTypes.tsv


# Step 4: Detection of somatic mutations
step4=Step4_VariantCalling
output_dir4=$output_dir/$step4
mkdir -p $output_dir4

docker run -v $output_dir:$mount_outdir -it scomatic BaseCellCalling.step1.py \
	--infile $mount_outdir/$step3/${sample}.BaseCellCounts.AllCellTypes.tsv \
	--outfile $mount_outdir/$step4/${sample} \
	--ref $REF

docker run -v $output_dir:$mount_outdir -it scomatic BaseCellCalling.step2.py \
	--infile $mount_outdir/$step4/${sample}.calling.step1.tsv \
	--outfile $mount_outdir/$step4/${sample} \
	--editing $editing \
	--pon $PON

docker run -v $output_dir:$mount_outdir -it scomatic bash -c "bedtools intersect -header -a $mount_outdir/$step4/${sample}.calling.step2.tsv -b $BOI | awk '\$1 ~ /^#/ || \$6 == \"PASS\"' > $mount_outdir/$step4/${sample}.calling.step2.pass.tsv"

# Computing the number of callable sites per cell type

step5=CellTypeCallableSites
output_dir5=$output_dir/$step5
mkdir -p $output_dir5

docker run -v $output_dir:$mount_outdir -it scomatic GetAllCallableSites.py \
	--infile $mount_outdir/$step4/${sample}.calling.step1.tsv  \
	--outfile $mount_outdir/$step5/${sample} \
	--max_cov 150 --min_cell_types 2

# Computing the number of callable sites per cell
step6=UniqueCellCallableSites
output_dir6=$output_dir/$step6
mkdir -p $output_dir6

for bam in $(ls -d $output_dir1/*bam);do  
	cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
	echo $cell_type
	
	fname=$(basename $bam)

	temp=$output_dir6/temp_${cell_type}
	mkdir -p $temp
	mount_temp=$mount_outdir/$step6/temp_${cell_type}

	docker run -v $output_dir:$mount_outdir -it scomatic SitesPerCell.py \
	  --bam $mount_outdir/$step1/$fname \
		--infile $mount_outdir/$step4/${sample}.calling.step1.tsv   \
		--ref $REF \
		--out_folder $mount_outdir/$step6 \
		--tmp_dir $mount_temp \
		--nprocs 1
done

# Computing the genotype for each cell at the variant sites

step7=SingleCellAlleles
output_dir7=$output_dir/$step7
mkdir -p $output_dir7

for bam in $(ls -d $output_dir1/*bam);do  
    cell_type=$(basename $bam | awk -F'.' '{print $(NF-1)}')
		fname=$(basename $bam)
    
    temp=$output_dir7/temp_${cell_type}
    mkdir -p $temp
		mount_temp=$mount_outdir/$step7/temp_${cell_type}

    docker run -v $output_dir:$mount_outdir -it scomatic SingleCellGenotype.py \
			--bam $mount_outdir/$step1/$fname \
			--infile ${mount_outdir}/$step4/${sample}.calling.step2.pass.tsv   \
			--nprocs 1   \
			--meta $META   \
			--outfile $mount_outdir/$step7/${cell_type}.single_cell_genotype.tsv  \
			--tmp_dir $mount_temp  \
			--ref $REF

    rm -rf $temp
done

# Computing the trinucleotide context background
step8=TrinucleotideContext
output_dir8=$output_dir/$step8
mkdir -p $output_dir8

docker run -v $output_dir:$mount_outdir -it scomatic bash -c "echo $mount_outdir/$step4/${sample}.calling.step1.tsv > $mount_outdir/$step8/step1_files.txt"

docker run -v $output_dir:$mount_outdir -it scomatic TrinucleotideContextBackground.py \
	--in_tsv $mount_outdir/$step8/step1_files.txt \
	--out_file $mount_outdir/$step8/TrinucleotideBackground.txt