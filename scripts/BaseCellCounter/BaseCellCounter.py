import pysam
import timeit
import pybedtools
import multiprocessing as mp
import argparse
import pandas as pd
import glob
import os
import sys
import time

def collect_result(result):
	if (result[1] != ''):
		VARIANTS = result[1]
		OUT = result[0]
		out = open(OUT,'w')
		out.write(VARIANTS)
		out.close()
		out_temp = result[0]


def concatenate_sort_temp_files_and_write(out_file, tmp_dir, ID):
	# Get the file paths
	all_files = glob.glob(tmp_dir + '/*.BaseCellCounts.temp')
	
	# Load as panda files
	if (len(all_files) > 0):
		
		## Organize files in dictionaries of chromosomes and starts
		Dictionary_of_files = {}
		for filename in all_files:
			basename = os.path.basename(filename)
			
			# Old way. Deprecated
			# basename = basename.split(".")
			# coordinates = basename[-3]
			
			# It can consider not standard chromosome nomenglatures (with dots)
			coordinates = basename.replace('.BaseCellCounts.temp','')

			CHROM, START, END = coordinates.split("__")
			
			START = int(START)
			
			if (CHROM not in Dictionary_of_files):
				Dictionary_of_files[CHROM] = {}
				Dictionary_of_files[CHROM][START] = filename
			else:
				Dictionary_of_files[CHROM][START] = filename
		

		## Write in the final output file
		out = open(out_file,'w')
		Header=['#CHROM','POS', 'REF', 'INFO', str(ID)]
		date = time.strftime("%d/%m/%Y")## dd/mm/yyyy format
		DATE="##fileDate=%s" % date
		CONCEPTS="""##INFO=DP,Description="Depth of coverage">\n##INFO=NC,Description="Number of different cells">\n##INFO=CC,Description="Cell counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BC,Description="Base counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BQ,Description="Base quality sums [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BCf,Description="Base counts in forward reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BCr,Description="Base counts in reverse reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">"""

		out.write(DATE + '\n')
		out.write(CONCEPTS + '\n')
		out.write('\t'.join(Header)+'\n')
		
		# Move thrugh filenanes by sorted coordinates
		for chrom in sorted(Dictionary_of_files.keys()):
			for start in sorted(Dictionary_of_files[chrom].keys()):
				filename = Dictionary_of_files[chrom][start]
				
				with open(filename, 'r') as f:
					out.write(f.read())
					out.write('\n')
		
				# Remove temp file
				os.remove(filename)
		out.close()
			
	
	else:
		# If temp files not found, print message
		print ('No temporary files found')

def MakeWindows(CONTIG, FASTA, bed, bed_out, window):
	## Pass 1. Get bed file and focus on the selected chromosome (if specified)
	if (bed == ''):
		# We get the bed file based on all coordenates from reference fasta file
		inFasta = pysam.FastaFile(FASTA)
		CONTIG_Names = inFasta.references
		LIST = [ (x, 1, inFasta.get_reference_length(x)) for x in CONTIG_Names]
		a = pybedtools.BedTool(LIST)
	else:
		# If bed file is provided, use this one
		a = pybedtools.BedTool(bed)
	
	# Focus only on the chromosome of interest (if provided)
	if (CONTIG != 'all'):
		a2 = a.filter(lambda b: b.chrom == CONTIG)
	else:
		a2 = a
	
	## Pass2. Intersect out (subtract non-desired regions) (if provided)
	if (bed_out != ''):
		b = pybedtools.BedTool(bed_out)
		a3 = a2.subtract(b)
	else:
		a3 = a2
	
	## Pass3. Makewindows in final bed file based on the window sizes specified as argument
	final_bed = a3.window_maker(a3,w=window)
	inFasta.close()
	return(final_bed)

def BaseCount(LIST, REF_BASE):
	Bases=['A','C','T','G','N']
	
	# Dictinary with base counts
	NUCLEOTIDES = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'I': 0, 'D' : 0, 'N': 0, 'R' : 0, 'Other' : 0}
	
	# Update dictionary counts
	for x in LIST:
		if x.upper() in Bases:
			NUCLEOTIDES[x.upper()] += 1
		elif '-' in x:
			NUCLEOTIDES['D'] += 1
		elif '+' in x:
			NUCLEOTIDES['I'] += 1
		else:
			NUCLEOTIDES['Other'] += 1
	
	# Calculate Alternative count of the most frequent alternative allele    
	Alternative = ['A','C','T','G','I', 'D']
	
	ALTERNATIVE = {k:v for k,v in NUCLEOTIDES.items() if k != REF_BASE.upper() and k in Alternative}

	AC = max(ALTERNATIVE.items(), key=lambda x: x[1])[1]
	
	if (AC > 0):
		listOfKeys = list()
		# Iterate over all the items in dictionary to find keys with max value
		for key, value in ALTERNATIVE.items():
			if value == AC:
				listOfKeys.append(key)
		
		MAX_ALT = ','.join(sorted(listOfKeys))
	else:
		MAX_ALT = '.'
		
	return ([NUCLEOTIDES,MAX_ALT,AC])

def EasyReadPileup(LIST, REF_BASE):
	Bases = set(['A','C','T','G','N'])
	
	AC = 0

	# Create a new list based on pileup reading
	NEW_LIST = []
	for x in LIST:
		LEN = len(x)
		UPPER = x.upper()
		if UPPER in Bases:
			NEW_LIST.append(UPPER)
			# Check for Alt counts
			if (UPPER != REF_BASE):
				AC = AC + 1
		elif LEN > 1 and x[1] == '-':
			D = 'D'
			NEW_LIST.append(D)
			AC = AC + 1
		elif LEN > 1 and x[1] == '+':
			I = 'I'
			NEW_LIST.append(I)
			AC = AC + 1
		elif x == '*':
			NEW_LIST.append('O')
		else:
			NEW_LIST.append('NA')
			
	return (NEW_LIST,AC)

def run_interval(interval, BAM, FASTA, MIN_COV, MIN_CC,MIN_AF, MIN_AC, tmp_dir, BQ, MQ):
	
	# Coordinates to analyse
	CHROM  =  interval[0]
	START = int(interval[1])
	END = int(interval[2])
	
	# Get pileup read counts from coordinates
	bam = pysam.AlignmentFile(BAM)
	i = bam.pileup(CHROM, START, END, min_base_quality = BQ, min_mapping_quality = MQ, ignore_overlaps = False)
	
	# Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
	inFasta = pysam.FastaFile(FASTA)

	# Run it for each position in pileup
	POSITIONS = []
	for p in i:
		POS=p.pos
		if POS >= START and POS < END:
			# Get reference base from fasta file
			ref_base = inFasta.fetch(CHROM, POS, POS+1)
			ref_base = ref_base.upper()
			
			# Get coverage
			DP = p.get_num_aligned()
			
			CELLS = []

			# Run only if coverage is more than minimum (arg parameter)
			if (DP >= MIN_COV and ref_base.upper() != 'N'):
				
				# Get pileup info
				READS = p.get_query_names()
				QUALITIES = p. get_query_qualities()
				PILEUP_LIST = p.get_query_sequences(mark_matches=True, add_indels=True)
				NEW_PILEUP_LIST,AC = EasyReadPileup(PILEUP_LIST, ref_base)

				# As intronic regions are also covered, we don't count them for the coverage filter
				CHECKING_CALLABLE_SITES = [x for x in NEW_PILEUP_LIST if x != 'NA']
				if (len(CHECKING_CALLABLE_SITES) < MIN_COV or AC < MIN_AC):
					continue

				# Get reads info
				reads = p.pileups

				# Dictionary counters
				BASE_COUNTS = {'A': 0, 'C': 0, 'T': 0,'G': 0, 'D': 0, 'I': 0, 'N': 0, 'O' : 0}
				BASE_QUALITIES = {'A': 0, 'C': 0, 'T': 0,'G': 0, 'D': 0, 'I': 0, 'N': 0, 'O' : 0}
				CELL_COUNTS = {'A': [], 'C': [], 'T': [],'G': [], 'D': [], 'I': [], 'N': [], 'O' : []}
				BASE_COUNTS_f = {'A': 0, 'C': 0, 'T': 0,'G': 0, 'D': 0, 'I': 0, 'N': 0, 'O' : 0}
				BASE_COUNTS_r = {'A': 0, 'C': 0, 'T': 0,'G': 0, 'D': 0, 'I': 0, 'N': 0, 'O' : 0}

				# Check this base in each read
				count = 0
				for read_i in range(0, len(READS)):
					
					read_alignment = reads[read_i].alignment

					try:
						barcode = read_alignment.opt("CB")
					except:
						continue

					# Fix barcode
					barcode = barcode.split("-")[0]

					# Remove low quality reads: Seconday alignments, duplicate reads, supplementary alignments
					if read_alignment.is_secondary == False and read_alignment.is_duplicate == False and read_alignment.is_supplementary == False:
						# Read info for the analysing position
						read = READS[read_i]
						
						# Get exact base
						base = NEW_PILEUP_LIST[read_i]
						base2 = PILEUP_LIST[read_i] 
						bq = QUALITIES[read_i]
						
						if base in BASE_COUNTS.keys(): # Usually, sites between exones (intronic sites) are marked as > or < . Therefore, we want to ignore them
							
							# Quality count
							count = count + 1

							# Update base count
							BASE_COUNTS[base] = BASE_COUNTS[base] + 1

							# Update base count
							BASE_QUALITIES[base] = BASE_QUALITIES[base] + bq

							# Strand info
							if read_alignment.is_reverse == True:
								BASE_COUNTS_r[base] = BASE_COUNTS_r[base] + 1
							else:
								BASE_COUNTS_f[base] = BASE_COUNTS_f[base] + 1

							# Update base count
							CELL_COUNTS[base].append(barcode)

							# Append cell barcode to check total number of cells
							CELLS.append(barcode)


				if (count >= MIN_COV):
					CELL_COUNTS2 = {x:len(set(CELL_COUNTS[x])) for x in CELL_COUNTS.keys()}


					# Pysam provides a 0-based coordinate. We must sum 1 to this value
					CHROM = str(CHROM)
					POS_print = str(POS + 1)
					REF = str(ref_base)

					DP = str(count)
					NC = len(set(CELLS))

					if (NC >= MIN_CC):
						NC = str(NC)

						info_field = '|'.join(['DP','NC','CC','BC','BQ','BCf','BCr'])

						# Lines to print
						Alleles = ['A','C','T','G','I', 'D']

						BASE_COUNTS_str = ':'.join([str(BASE_COUNTS[x]) for x in  Alleles])
						BASE_QUALITIES_str = ':'.join([str(BASE_QUALITIES[x]) for x in  Alleles])
						BASE_COUNTS_f_str = ':'.join([str(BASE_COUNTS_f[x]) for x in  Alleles])
						BASE_COUNTS_r_str = ':'.join([str(BASE_COUNTS_r[x]) for x in  Alleles])
						CELL_COUNTS_str = ':'.join([str(CELL_COUNTS2[x]) for x in  Alleles])

						INFO = "|".join([DP,NC,CELL_COUNTS_str,BASE_COUNTS_str,BASE_QUALITIES_str,BASE_COUNTS_f_str,BASE_COUNTS_r_str])
						LINE = '\t'.join([CHROM, POS_print, REF,info_field,INFO])

						# Save results in list of positions
						POSITIONS.append(LINE)
			
	inFasta.close()
	bam.close()
	
	# Return list of positions
	ID = '__'.join([str(CHROM), str(START), str(END)])
	out_temp = tmp_dir + '/' + ID + '.BaseCellCounts.temp'
	return([out_temp,'\n'.join(POSITIONS)])


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to obtain a list of base and cell counts in scRNA bam file')
	parser.add_argument('--bam', type=str, default=1, help='BAM file to be analysed', required = True)
	parser.add_argument('--ref', type=str, default=1, help='Path to reference genome version. *fai must be available in the same directory as the reference genome file', required = True)
	parser.add_argument('--chrom', type=str, help='Chromosome to be analysed. --chrom all to analyse all chromosomes', required = True)
	parser.add_argument('--out_folder', default = '.', help='Out folder', required = False)
	parser.add_argument('--id', help='Prefix used to name output file. If provided, please conform with the following format: *.[cell_type] . Example: sample1.t_cell. If not provided, the basename of the BAM file will be used.', required = False)
	parser.add_argument('--nprocs',default = 1, help='Number of processes [Default: 1]',required=False,type = int)
	parser.add_argument('--bin', type=int, default=50000, help='Bin size for running the analysis [Default: 50000]', required = False)
	parser.add_argument('--bed', type=str, default='', help='Regions to focus the analysis on. Three-column bed file listing the chromosome, start and end for those regions to be analysed.', required = False)
	parser.add_argument('--bed_out', type=str, default='', help='Regions to ignore in the analysis. Three-column bed file listing the chromosome, start and end for those regions to be ignored.', required = False)
	parser.add_argument('--min_ac', type=int, default = 0, help='Minimum number of reads supporting the alternative allele required to consider a genomic site for mutation calling. Default: 0', required = False)
	parser.add_argument('--min_af', type=float, default = 0, help='Minimum alternative allele fraction required to consider a genomic site for mutation calling. Default = 0', required = False)
	parser.add_argument('--min_dp', type=int, default = 5, help='Minimum depth of coverage required to consider a genomic site for mutation calling. Default: 5', required = False)
	parser.add_argument('--min_cc', type=int, default = 5, help='Minimum number of cells required to consider a genomic site for mutation calling. Default: 5', required = False)
	parser.add_argument('--min_bq', type=int, default = 20, help='Minimum base quality to compute allele counts. Default: 20', required = False)
	parser.add_argument('--min_mq', type=int, default = 255, help='Minimum mapping quality required to consider a read for analysis. Default: 255', required = False)
	parser.add_argument('--tmp_dir', type=str, default = '.', help='Path to a directory to be used to store temporary files during processing', required = False)

	return (parser)

def main():


	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	BAM = args.bam
	FASTA = args.ref
	CORE = args.nprocs
	OUTF = args.out_folder
	ID = args.id
	CONTIG = args.chrom
	BIN = args.bin
	bed = args.bed
	bed_out = args.bed_out
	MIN_AC = args.min_ac
	MIN_AF = args.min_af
	MIN_COV = args.min_dp
	MIN_CC = args.min_cc
	MIN_BQ = args.min_bq
	MIN_MQ = args.min_mq
	tmp_dir = args.tmp_dir

	# 0. If out file is empty, take the name of the bam file
	if (ID == None):
		ID = os.path.basename(BAM)
		ID = ID.replace(".bam","")

	# Set outfile name
	out_file = OUTF + '/' + str(ID) + ".tsv"
	print("Outfile: " , out_file ,  "\n") 

	# 1. Check if temp dir is created
	if (tmp_dir != '.'):
		try:
			# Create target Directory
			os.mkdir(tmp_dir)
			print("Directory " , tmp_dir ,  " created\n") 
		except FileExistsError:
			print("Directory " , tmp_dir ,  " already exists\n")
	else:
		print("Not temp directory specified, using working directory as temp") 

	# 2. Create bed file and windows
	BED = MakeWindows(CONTIG, FASTA, bed, bed_out, BIN)
	
	# 3. Code to run in parallel all bins
	if (CORE > 1):
		pool = mp.Pool(CORE)
		
		# Step 3.1: Use loop to parallelize
		for row in BED:
			# This funtion writes in temp files the results
			pool.apply_async(run_interval, args=(row, BAM, FASTA, MIN_COV, MIN_CC, MIN_AF, MIN_AC, tmp_dir, MIN_BQ, MIN_MQ), callback=collect_result)
				   
		# Step 3.2: Close Pool and let all the processes complete    
		pool.close()
		pool.join()
	else:
		for row in BED:
			# This funtion writes in temp files the results
			collect_result(run_interval(row, BAM, FASTA, MIN_COV, MIN_CC, MIN_AF, MIN_AC,tmp_dir, MIN_BQ, MIN_MQ))
	
	# 4. Write final file
	concatenate_sort_temp_files_and_write(out_file, tmp_dir, ID)


if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')


