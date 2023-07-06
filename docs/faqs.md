# FAQs - Frequently asked questions

## 1. Are the SComatic parameters for scATAC-seq data the same as for scRNA-seq data?
No, they are different. ATAC-seq data is a DNA-based approach. Therefore, there are differences in the way of processing the bam files. These are the parameters that you should change in comparison to the scRNA-seq workflow:
- **Step 1**: remove the _--max_nM_ and _--max_NH_ parameters, set the mapping quality filter to _--min_MQ_ 30 .
- **Step 2**: set the mapping quality filter to _--min_mq_ 30
- **Step 4.2**: Remove the _--editing_ parameter, and _--pon_ altered to point at the scATACseq PoN provided in this GitHub repo (or custom one)

## 2. How can we perform the variant annotation with the SComatic output?
In our manuscript, we have annotated all our variants using [annovar](https://annovar.openbioinformatics.org/en/latest/). With a couple of command lines, we can adapt our variant calling output for this annotation tool. Using our example data, we should run these lines:

**1.  Preparing SComatic output for annovar**

```
ANNOVAR=/path/to/annovar
hummandb=/path/to/annovar_humandb_hg38

grep -v '#' Example.calling.step2.pass.tsv |  tr '\t' '-' | awk -F'-' -v OFS='\t' '{print $1,$2,$3,$4,$5,$0}' > sample.variants.avinput
```

**2. Annotate variants using the annovar-ready variant file**

```
perl $ANNOVAR/table_annovar.pl sample.variants.avinput \
                $hummandb -buildver hg38 \
                -out sample.variants.annovar \
                -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a,gnomad_genome -operation g,r,f,f,f,f \
                -nastring . -csvout -polish --otherinfo
```

This step will generate a comma-separated file (.csv) called _sample.variants.annovar.hg38_multianno.csv_ , which can be easily opened and processed in your more desired software (R, excel...). The output should look like this:

```
Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,cytoBand,ExAC_ALL,ExAC_AFR,ExAC_AMR,ExAC_EAS,ExAC_FIN,ExAC_NFE,ExAC_OTH,ExAC_SAS,avsnp147,SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_score,MutationAssessor_pred,FATHMM_score,FATHMM_pred,PROVEAN_score,PROVEAN_pred,VEST3_score,CADD_raw,CADD_phred,DANN_score,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,integrated_fitCons_score,integrated_confidence_value,GERP++_RS,phyloP7way_vertebrate,phyloP20way_mammalian,phastCons7way_vertebrate,phastCons20way_mammalian,SiPhy_29way_logOdds,gnomAD_genome_ALL,gnomAD_genome_AFR,gnomAD_genome_AMR,gnomAD_genome_ASJ,gnomAD_genome_EAS,gnomAD_genome_FIN,gnomAD_genome_NFE,gnomAD_genome_OTH,Otherinfo
chr10,29559501,29559501,A,T,"intronic","SVIL",.,.,.,"10p11.23",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-29559501-29559501-A-T-PASS-Epithelial_cells-ATTCA-TTGAT-1-12-5-4-2-0.3333-0.4-0.0-0.0001-2-2-0;38;1-0;16;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-12|5|3:0:2:0:0:0|8:0:4:0:0:0|306:0:164:0:0:0|0:0:0:0:0:0|8:0:4:0:0:0-30|13|13:0:0:0:0:0|30:0:0:0:0:0|1189:0:0:0:0:0|0:0:0:0:0:0|30:0:0:0:0:0"
chr10,73413109,73413109,G,A,"intronic","ANXA7",.,.,.,"10q22.2",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-73413109-73413109-G-A-PASS-Epithelial_cells-TGGAG-CATGG-1-12-8-4-3-0.3333-0.375-0.0-0.0-2-2-0;15;1-0;11;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-12|8|3:0:0:5:0:0|4:0:0:8:0:0|160:0:0:302:0:0|0:0:0:0:0:0|4:0:0:8:0:0-7|6|0:0:0:6:0:0|0:0:0:7:0:0|0:0:0:270:0:0|0:0:0:0:0:0|0:0:0:7:0:0"
chr10,89413978,89413978,G,A,"upstream","IFIT5","dist=590",.,.,"10q23.31",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-89413978-89413978-G-A-PASS-Epithelial_cells-GTATG-GTGTT-1-18-11-7-4-0.3889-0.3636-0.0-0.0-2-2-0;20;1-0;12;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-18|11|4:0:0:7:0:0|7:0:0:11:0:0|270:0:0:417:0:0|0:0:0:0:0:0|7:0:0:11:0:0-9|5|0:0:0:5:0:0|0:0:0:9:0:0|0:0:0:339:0:0|0:0:0:0:0:0|0:0:0:9:0:0"
chr10,97332670,97332670,G,A,"UTR3","FRAT2","NM_012083:c.*1201C>T",.,.,"10q24.1",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-97332670-97332670-G-A-PASS-Epithelial_cells-TTGGG-GTGGC-1-12-11-6-6-0.5-0.5455-0.0-0.0-2-2-0;12;1-0;11;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-6|6|0:0:0:6:0:0|0:0:0:6:0:0|0:0:0:237:0:0|0:0:0:0:0:0|0:0:0:6:0:0-12|11|6:0:0:5:0:0|6:0:0:6:0:0|229:0:0:225:0:0|0:0:0:0:0:0|6:0:0:6:0:0-NA"
chr10,126926205,126926205,G,A,"intronic","DOCK1",.,.,.,"10q26.2",.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,.,"chr10-126926205-126926205-G-A-PASS-Epithelial_cells-GATCC-GTGAC-1-35-17-15-6-0.4286-0.3529-0.0-0.0-2-2-0;35;1-0;18;1-.-PASS-DP|NC|CC|BC|BQ|BCf|BCr-NA-35|17|6:0:0:12:0:0|15:0:0:20:0:0|603:0:0:786:0:0|15:0:0:20:0:0|0:0:0:0:0:0-15|7|0:0:0:7:0:0|0:0:0:15:0:0|0:0:0:590:0:0|0:0:0:15:0:0|0:0:0:0:0:0"
```

As you will notice, the first few columns of the file have the annotation information provided by _annovar_ (Gene, region, impact, Gnomad allele frequencies...). Importantly, the final column of the csv file (_Otherinfo_) shows all the info found in the original _Example.calling.step2.pass.tsv_, but separated by the symbol "-" . The column names of each one of these _Otherinfo_ items are the same as the ones found in the _Example.calling.step2.pass.tsv_ : 

```
> grep '^#CHROM' Example.calling.step2.pass.tsv

#CHROM	Start	End	REF	ALT	FILTER	Cell_types	Up_context	Down_context	N_ALT	Dp	Nc	Bc	Cc	VAF	CCF	BCp	CCp	Cell_types_min_BC	Cell_types_min_CC	Rest_BC	Rest_CC	Fisher_p	Cell_type_Filter	INFO	Myeloids	Epithelial_cells	Stromal_cells

```

Of course, depending on the columns you want to keep for the annotation process, you might slightly change the command line in the _step 1_ described above, as well as the parameters provided in the annovar computation. The parameters provided here are the ones used in our manuscript. 

## 3. Can SComatic work with other types of PoN files?
Yes, SComatic can deal with PoNs obtained by other tools.

The *--pon* parameter permits working with different types of formats/files depending on the availability and quantity of non-neoplastic samples. These are the main options: 

* Using the [PoNs](/PoNs/) provided in this repository, which were computed using the data described in the main manuscript of SComatic.
* Using your custom PoNs, which can be built using [SComatic modules](/docs/pon.md).
* TSV file listing the mutations detected by running SComatic (output of *Detection of somatic mutations: Step 4.2*) on scRNA-seq data from e.g., matched non-neoplastic cells.
* Custom PoN in VCF format (unzipped) generated by running a variant caller, such as *GATK-HaplotypeCaller*, on DNA sequencing (WES or WGS) data, such as a matched normal sample or a set of unmatched germline samples.
* Using a PoN in VCF format (unzipped) generated by other tools like [*GATK*](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-)





