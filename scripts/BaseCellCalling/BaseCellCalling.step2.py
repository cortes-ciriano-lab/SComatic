#!/usr/bin/python

import numpy as np
import timeit
import os
import math
import pysam
import argparse
from scipy.stats import betabinom
import scipy.stats as stats
import pandas as pd
import subprocess



def variant_calling_step2(file,distance,editing,pon,window,outfile):

	#---
	# Command to focus only on candidate sites
	#---
	# Temporary file to work with
	file_temp = file + '.temp'
	command = "awk -F'\\t' '{if (($1 ~ /^#/) || ($6 != \".\")) {print $0}}' %s > %s " % (file,file_temp)

	# Submit linux command
	try:
		subprocess.run(command, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)

	#---
	# Run extra filters
	#---

	# Build RNA-diting and PoN dictionaries
	EDITING_DICT = build_dict(editing,window)
	PON_DICT = build_dict(pon,window)	

	# Run extra filters
	outfile= open(outfile,'w')
	current_chr = 0
	LIST = list()
	with open(file_temp, 'r') as f:
		for line in f:
			if line.startswith('#'):
				outfile.write(line)
			else:
				line = line.rstrip('\n')
				elements = line.split('\t')
				
				# We append new candidate sites in a list of max three candidates to check the distance with the variant up- and down-stream in the vcf-like file
				LIST.append(elements)

				# As soon as we have 3 sites in the list, we proceed with the filtering
				if len(LIST) == 3:
					
					# Don't forget first candidate of the file
					if current_chr == 0:
						candidate = LIST[0]
						current_chr = candidate[0] # First chromosome in file. Flag to start the analysis

						# Run extra filters for each candidate site
						NEW_LINE = GetExtraFilters(LIST,candidate,distance,EDITING_DICT,PON_DICT,window)
						outfile.write(NEW_LINE)

					# Candidate in the middle of candidate context (in a list of 3)
					candidate = LIST[1]
					NEW_LINE = GetExtraFilters(LIST,candidate,distance,EDITING_DICT,PON_DICT,window)
					outfile.write(NEW_LINE)

					# Remove first element
					LIST.pop(0)

		# Avoid losing sites when there are less than 3 candidate sites
		if current_chr == 0:
			for candidate in LIST:
				NEW_LINE = GetExtraFilters(LIST,candidate,distance,EDITING_DICT,PON_DICT,window)
				outfile.write(NEW_LINE)

		# Avoid losing last sites from the file
		elif (len(LIST) > 1):
			candidate = LIST[1]
			NEW_LINE = GetExtraFilters(LIST,candidate,distance,EDITING_DICT,PON_DICT,window)
			outfile.write(NEW_LINE)
	outfile.close()

	# Remove temp file
	os.remove(file_temp)


def GetExtraFilters(LIST,candidate,distance,EDITING_DICT,PON_DICT,window):
	# Get info from candidate variant
	candidate_chr = candidate[0]
	candidate_pos = int(candidate[1])
	
	current_chr = candidate_chr

	# Filter column
	FILTER = candidate[5]
	
	# Edit based on distance between variants, listed in editing site or listed in PoN list
	if (FILTER != "."):

		# Get if there are variants close to our candidarte variant
		# Return how many of the neighbours are too close	
		CLOSE = len([x[1] for x in LIST if x[0] == candidate_chr and int(x[1]) != candidate_pos and abs(int(x[1])-candidate_pos) <= distance])

		# Check for editing
		WIND = math.floor(candidate_pos / float(window))
		try: 
			EDITING = candidate_pos in EDITING_DICT[current_chr ][WIND]
		except:
			EDITING = False

		# Check for PoN
		try: 
			PON = candidate_pos in PON_DICT[current_chr ][WIND]
		except:
			PON = False


		# Check if there are potential variants close to the candidate site 	
		if CLOSE > 0 or EDITING == True or PON == True:
			# Editing filter
			if EDITING == True:
				if (FILTER == 'PASS'):
					FILTER  = 'RNA_editing_db'
				else:
					FILTER = FILTER + ',RNA_editing_db'

			# Close variants filter
			if (CLOSE > 0):
				if (FILTER == 'PASS'):
					FILTER  = 'Clustered'
				else:
					FILTER = FILTER + ',Clustered'

			# Editing filter
			if PON == True:
				if (FILTER == 'PASS'):
					FILTER  = 'PoN'
				else:
					FILTER = FILTER + ',PoN'

			candidate[5] = FILTER
	
	# Prepare line for printing
	candidate = '\t'.join(candidate) + '\n'
	return(candidate)

def build_dict(editing,window):
	DICT_editing = {}
	try:
		with open(editing, 'r') as f:
			for line in f:
				if not line.startswith('#'):
					elements = line.split('\t')

					# Coordinates
					CHROM = elements[0]
					POS = int(elements[1])

					WIND = math.floor(POS / float(window))

					if not CHROM in DICT_editing.keys():
						DICT_editing[CHROM] = {}
						DICT_editing[CHROM][WIND] = set([POS])
					else:
						if not WIND in DICT_editing[CHROM].keys():
							DICT_editing[CHROM][WIND] = set([POS])
						else:
							DICT_editing[CHROM][WIND].update([POS])
	except:
		DICT_editing = {}
	return (DICT_editing)


def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to perform the scRNA somatic variant calling')
	parser.add_argument('--infile', type=str, help='Input file with all samples merged in a single tsv', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--editing', type=str, help='RNA editing file to be used to remove RNA-diting sites', required = False)
	parser.add_argument('--pon', type=str, help='Panel of normals (PoN) file to be used to remove germline polymorphisms and recurrent artefacts', required = False)	
	parser.add_argument('--min_distance', type=int, default = 5, help='Minimum distance allowed between potential somatic variants [Default: 5]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	editing = args.editing
	pon = args.pon
	distance = args.min_distance

	# How to split genome (window sizes)
	window = 20000

	# 1. Variant calling
	print ('\n------------------------------')
	print ('Variant calling')
	print ('------------------------------\n')

	# 1.2: Step 2: Add distance, editing and PoN filters
	print ('\n- Variant calling step 2\n')
	print("	> Editing file used: " , editing)
	print("	> PoN file used: " , pon)
	outfile2 = outfile + '.calling.step2.tsv'
	variant_calling_step2(infile,distance,editing,pon,window,outfile2)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')


