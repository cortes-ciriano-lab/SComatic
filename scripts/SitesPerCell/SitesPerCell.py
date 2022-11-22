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
import subprocess
import shutil
import numpy as np
from collections import Counter

def collect_result(result):
	if (result[1] != ''):
		VARIANTS = result[1]
		OUT = result[0]
		np.save(OUT,VARIANTS)

def concatenate_sort_temp_files_and_write(out_file, tmp_dir, ID):
	# Get the file paths
	all_files = glob.glob(tmp_dir + '/*.SitesPerCell.npy')
	pattern=tmp_dir + '/*.SitesPerCell.npy'
	
	final_dict = dict()
	for npy in all_files:
		a = np.load(npy,allow_pickle=True)
		b = a.item()

		# Collapse and sum all cell counts
		for cb in b.keys():
			c = b[cb]
			final_dict[cb] = final_dict.get(cb, 0) + c

		# Remove file
		os.remove(npy)

	# Transform to data frame	
	df = pd.DataFrame.from_dict(final_dict, orient='index')
	df.index.name = 'CB'
	df.reset_index(inplace=True)
	df.rename({0: 'SitesPerCell'}, axis=1, inplace=True)
	df.to_csv(out_file, index=False)


def MakeWindows(bed, window):
	## Makewindows in bed file based on the window sizes specified as argument
	a = pybedtools.BedTool(bed)
	final_bed = a.window_maker(a,w=window)

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

def run_interval(interval,sites,BAM, FASTA, MIN_COV, MIN_CC, tmp_dir, BQ, MQ):
	
	# Coordinates to analyse
	CHROM,START,END = interval.split("_")
	START = int(START)-1
	END = int(END)+1
	
	# Get pileup read counts from coordinates
	bam = pysam.AlignmentFile(BAM)
	i = bam.pileup(CHROM, START, END, min_base_quality = BQ, min_mapping_quality = MQ, ignore_overlaps = False)
	
	# Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
	inFasta = pysam.FastaFile(FASTA)

	# Run it for each position in pileup
	CELLS_all = []
	for p in i:
		POS=p.pos

		if POS in sites:
			# Get reference base from fasta file
			ref_base = inFasta.fetch(CHROM, POS, POS+1)
			
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

				# Get reads info
				reads = p.pileups

				# Dictionary counters
				BASE_COUNTS = ['A', 'C', 'T','G', 'D', 'I']

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

						# Get exact base
						base = NEW_PILEUP_LIST[read_i]
						
						if base in BASE_COUNTS: # Usually, sites between exones (intronic sites) are marked as > or < . Therefore, we want to ignore them
						
							# Append cell barcode to check total number of cells
							CELLS.append(barcode)

						count=count+1

				if (count >= MIN_COV):
					CELLS = list(set(CELLS))


					# Pysam provides a 0-based coordinate. We must sum 1 to this value
					CHROM = str(CHROM)
					POS_print = str(POS + 1)
					REF = str(ref_base)

					DP = str(count)
					NC = len(set(CELLS))

					if (len(CELLS) >= MIN_CC):
						# Save results in list of positions
						CELLS_all.extend(CELLS)
			
	inFasta.close()
	bam.close()
	
	counts = dict()
	for i in CELLS_all:
		counts[i] = counts.get(i, 0) + 1

	# Return list of positions
	ID = '_'.join([str(CHROM), str(START), str(END)])
	out_temp = tmp_dir + '/' + ID + '.SitesPerCell.npy'
	return([out_temp,counts])


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to calculate the number of callable sites per unique cell')
	parser.add_argument('--bam', type=str, default=1, help='Tumor bam file to be analysed', required = True)
	parser.add_argument('--ref', type=str, default=1, help='Reference genome. *fai must be available in the same folder as reference', required = True)
	parser.add_argument('--infile', type=str, default='', help='Base calling file (obtained by BaseCellCalling.step1.py)', required = False)
	parser.add_argument('--min_ct1', type=int, default = 2, help='Minimum number of cell types with enough reads to consider a genomic site. Default = 2', required = False)
	parser.add_argument('--min_ct2', type=int, default = 2, help='Minimum number of cell types with enough unique cell counts to consider a genomic site. Default = 2', required = False)
	parser.add_argument('--out_folder', default = '.', help='Out folder', required = False)
	parser.add_argument('--id', help='Prefix of out file. If provided, please use next format: *.[cell_type] . Example: sample1.t_cell. If not provided, it is taken from bam file', required = False)
	parser.add_argument('--nprocs',default = 1, help='Number of processes [Default = 1]',required=False,type = int)
	parser.add_argument('--bin', type=int, default=50000, help='Bin size for running the analysis [Default 50000]', required = False)
	parser.add_argument('--min_dp', type=int, default = 5, help='Minimum coverage to consider the genomic site. Default = 5', required = False)
	parser.add_argument('--min_cc', type=int, default = 5, help='Minimum number of cells required to consider a genomic site. Default = 5', required = False)
	parser.add_argument('--min_bq', type=int, default = 20, help='Minimum base quality permited for the base counts. Default = 20', required = False)
	parser.add_argument('--min_mq', type=int, default = 255, help='Minimum mapping quality required to analyse read. Default = 255', required = False)
	parser.add_argument('--tmp_dir', type=str, default = '.', help='Temporary folder for tmp files', required = False)

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
	BIN = args.bin
	in_tsv = args.infile
	MIN_CT1 = args.min_ct1
	MIN_CT2 = args.min_ct2
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
	out_file = OUTF + '/' + str(ID) + ".SitesPerCell.tsv"
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



	# Temp file and folder
	temp = tmp_dir + '/target.bed'

	# 1. Concatenate and create temp bed file
	print ('-----------------------------------------------------------')
	print ('1. Preparing temp file...')
	print ('-----------------------------------------------------------\n')

	command = "grep -v '^#' %s | grep -v 'chrM' | awk -F'\\t' -v OFS='\\t' '{if ($19 >= %s && $20 >= %s) {print $1,$3}}'| sort -T %s -k1,1 -k2,2n > %s" % (in_tsv,MIN_CT1, MIN_CT2,tmp_dir,temp)

	# Submit linux command
	try:
		subprocess.run(command, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)

	# Split the in_tsv file in windows
	DICT_sites = {}
	cur_chr = ''
	cur_pos = 0
	count = 0
	with open(temp,'r') as tmp:
		for line in tmp:
			line = line.rstrip('\n')
			CHROM, POS = line.split("\t")
			POS = int(POS)

			if (cur_chr == ''):	
				cur_chr = CHROM
				cur_pos = POS
				LIST=[]

			DIFF = POS-cur_pos

			if CHROM == cur_chr and DIFF < 10000 and count < 1000:
				LIST.append(POS-1) # Necessary to subtract one for pysam coordinates
				count = count + 1
			else:
				ID=str(cur_chr) + '_' + str(cur_pos) + '_' + str(POS)
				DICT_sites[ID] = set(LIST)

				# Re-start counting
				cur_chr = CHROM
				cur_pos = POS
				count=0
				LIST= [POS-1] # Necessary to subtract one for pysam coordinates

		# Save results with last iteration in file
		if (count > 0):
			ID=str(cur_chr) + '_' + str(cur_pos) + '_' + str(POS)
			DICT_sites[ID] = set(LIST)




	# 2. Checking bam file
	print ('-----------------------------------------------------------')
	print ('2. Checking bam file...')
	print ('-----------------------------------------------------------\n')

	if (CORE > 1):
		pool = mp.Pool(CORE)
		
		# Step 3.1: Use loop to parallelize
		for row in DICT_sites.keys():
			# This funtion writes in temp files the results
			pool.apply_async(run_interval, args=(row, DICT_sites[row], BAM, FASTA, MIN_COV, MIN_CC, tmp_dir, MIN_BQ, MIN_MQ), callback=collect_result)

		# Step 3.2: Close Pool and let all the processes complete    
		pool.close()
		pool.join()
	else:
		for row in DICT_sites.keys():
			# This funtion writes in temp files the results
			collect_result(run_interval(row, DICT_sites[row], BAM, FASTA, MIN_COV, MIN_CC, tmp_dir, MIN_BQ, MIN_MQ))
	

	# 4. Write final file
	concatenate_sort_temp_files_and_write(out_file, tmp_dir, ID)

	# Remove temp directory
	shutil.rmtree(tmp_dir, ignore_errors=True)

if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ("Computation time: " + str(round(stop - start)) + ' seconds')


