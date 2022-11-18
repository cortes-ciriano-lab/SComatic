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

def reinit_template(diff_cell_types,min_reads,min_cells):
	# Creates empty dicctionary for coverage counts
	TEMPLATE = {min_reads:0, min_cells:0, 5: 0, 10: 0, 20: 0,30: 0}
	INFORMATIVE_POSITIONS_TEMPLATE = {x:{'NC':TEMPLATE.copy(),'DP':TEMPLATE.copy()} for x in diff_cell_types}
	return(INFORMATIVE_POSITIONS_TEMPLATE, TEMPLATE)

def	callable_sites(file,min_cell_types,MAX):
	Alleles = ["A","C","T","G","I","D","N","O"]
	VARIANT_POSITIONS = {}
	CHR_INFORMATIVE_POSITIONS = {}

	cur_chr = 'z'

	counter = 0
	with open(file, 'r') as f:
		for line in f:
			if line.startswith('##'):
				pass

			else:
				if line.startswith('#CHROM') and counter == 0:
					line = line.rstrip('\n')

					# Create dictionary with cell types and indexes for each one
					elements = line.split('\t')
					cell_types_idx = {x : elements[x] for x in range(len(elements)) if x > 24} # First 4 columns are common across all lines ['#CHROM', 'Start','End', 'REF', 'INFO']. Up to the 25th, they are lines that are created during the first step of the variant calling

					diff_cell_types = list(cell_types_idx.values())

					for cell_type in diff_cell_types:
						CHR_INFORMATIVE_POSITIONS[cell_type] = {}

				else:
					counter = counter + 1
					line = line.rstrip('\n')
					
					# Create dictionary with cell types and indexes for each one
					elements = line.split('\t')

					# Chromosome
					CHROM = str(elements[0])
					if (CHROM != cur_chr):
						for cell_type in CHR_INFORMATIVE_POSITIONS.keys(): # Create templates (setting values to 0 for the counting)
							CHR_INFORMATIVE_POSITIONS[cell_type][CHROM] = {}
							CHR_INFORMATIVE_POSITIONS[cell_type][CHROM] = {x:{'DP':0,'NC':0} for x in range(0,MAX+1)}
						cur_chr = CHROM

					# Position
					POS = int(elements[1])

					# Reference allele
					REF = elements[3]

					# Number of cell types with enough coverage
					Cell_types_min_DP = int(elements[18])  
					# Number of cell types with enough coverage
					Cell_types_min_NC = int(elements[19])   


					# Filter by the minimum number of cell types with enough coverage and number of cells
					if (Cell_types_min_DP >= min_cell_types and Cell_types_min_NC >= min_cell_types):

						### Site variant calling level info
						FILTER = []

						# Computation for each cell type
						for cell_type_i in cell_types_idx.keys():
							cell_type = cell_types_idx[cell_type_i]
							INFO_i = elements[cell_type_i]


							if not INFO_i.startswith('NA'):
								DP,NC,CC,BC,BQ,BCf,BCr = INFO_i.split('|')
								
								# Set a maximum to be stored in the callable table (dictionary so far)
								if (int(DP) > MAX):
									DP = MAX
								else:
									DP = int(DP)
								if (int(NC) > MAX):
									NC = MAX
								else:
									NC = int(NC)

								# Add values to dictionary
								CHR_INFORMATIVE_POSITIONS[cell_type][CHROM][DP]['DP'] += 1 
								CHR_INFORMATIVE_POSITIONS[cell_type][CHROM][NC]['NC'] += 1

	return(CHR_INFORMATIVE_POSITIONS)



def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to calculate the number of callable sites per cell type')
	parser.add_argument('--infile', type=str, help='Tsv file obtained by BaseCellCalling.step1.py to be analysed', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--max_cov', type=int, default = 150, help='Maximum coverage to record in the callable sites table. Greater values will be collapsed to the provided one. [Default: 150]', required = False)
	parser.add_argument('--min_cell_types', type=int, default = 2, help='Minimum number of cell types with enough coverage/cells at a site to be considered as a callable [Default: 2]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	max_cov = args.max_cov
	min_cell_types = args.min_cell_types

	# 1. Variant calling
	print ('\n------------------------------')
	print ('Getting callable sites')
	print ('------------------------------\n')

	# Getting callable sites
	start = timeit.default_timer()
	informative_positions = callable_sites(infile,min_cell_types,max_cov)

	# Coverage report
	outfile_report_all = outfile + '.coverage_cell_count.report.tsv'
	outfile_report_all_chrom = outfile + '.coverage_cell_count.per_chromosome.report.tsv'
	

	# Coverage per chromosome
	out1 = open(outfile_report_all_chrom,'w')
	Header1=['Cell_types','CHROM','Cov', 'DP','NC']
	out1.write('\t'.join(Header1)+'\n')
	for cell_type in informative_positions.keys():
		for chrom in informative_positions[cell_type].keys():
			for COV in informative_positions[cell_type][chrom]:
				DP = informative_positions[cell_type][chrom][COV]['DP']
				NC = informative_positions[cell_type][chrom][COV]['NC']
				LINE = [str(cell_type),str(chrom),str(COV),str(DP),str(NC)]
				LINE = '\t'.join(LINE)
				out1.write(LINE+'\n')
	out1.close()


	# Coverage all chromosomes collapsed
	DF = pd.read_csv(outfile_report_all_chrom, delimiter = '\t')
	DFtotal = DF.groupby(['Cell_types','Cov']).agg(DP = ('DP',sum), NC = ('NC',sum)).reset_index()
	DFtotal.to_csv(outfile_report_all, index=False, sep = '\t')


#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')
