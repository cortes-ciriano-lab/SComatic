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

def	callable_sites(file,min_cells,min_reads,min_cell_types):
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
					cell_types_idx = {x : elements[x] for x in range(len(elements)) if x > 24} # First 4 columns are common across all lines ['#CHROM', 'Start','End', 'REF', 'INFO']. The rest are cell types

					diff_cell_types = list(cell_types_idx.values())

					# # Let's save the number of informative sites
					# TEMPLATE = {min_reads:0, min_cells:0, 5: 0, 10: 0, 20: 0,30: 0}
					# INFORMATIVE_POSITIONS_TEMPLATE = {elements[x]:{'NC':TEMPLATE.copy(),'DP':TEMPLATE.copy()} for x in range(len(elements)) if x > 3}

				else:
					counter = counter + 1
					line = line.rstrip('\n')
					
					# Create dictionary with cell types and indexes for each one
					elements = line.split('\t')

					# Chromosome
					CHROM = str(elements[0])
					if (CHROM != cur_chr):
						CHR_INFORMATIVE_POSITIONS[CHROM], TEMPLATE = reinit_template(diff_cell_types,min_reads,min_cells)
						cur_chr = CHROM

					# Position
					POS = int(elements[1])

					# Reference allele
					REF = elements[3]

					### Single cell type level info
					# All DPs
					All_DPs = {}
					# All CCs
					All_NCs = {}
					# Number of cell types with enough coverage
					Cell_types_min_DP = 0  
					# Number of cell types with enough coverage
					Cell_types_min_NC = 0  

					### Site variant calling level info
					FILTER = []

					# Computation for each cell type
					for cell_type_i in cell_types_idx.keys():
						cell_type = cell_types_idx[cell_type_i]
						INFO_i = elements[cell_type_i]

						if INFO_i.startswith('NA'):
							All_DPs[cell_type] = 0
							All_NCs[cell_type] = 0
 
						elif not INFO_i.startswith('NA'):
							DP,NC,CC,BC,BQ,BCf,BCr = INFO_i.split('|')
							
							All_DPs[cell_type] = int(DP)
							All_NCs[cell_type] = int(NC)

							if (int(DP) >= min_cells and int(NC) >= min_reads):
								# Number of cell types with min_BC
								Cell_types_min_DP = Cell_types_min_DP + 1
								# Number of cell types with min_CC
								Cell_types_min_NC = Cell_types_min_NC + 1

					# Save coverage and cell count information
					if (Cell_types_min_NC >= min_cell_types):
						# DPs
						for cell_type in All_DPs.keys():
							DP = All_DPs[cell_type]
							NC = All_NCs[cell_type]

							for x in TEMPLATE.keys():
								if DP >= x:
									CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['DP'][x] = CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['DP'][x] + 1
								if NC >= x: 
									CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['NC'][x] = CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['NC'][x] + 1

					if (counter % 1000000 == 0):
						print ('Getting callable sites per cell type: ' + str(counter) + ' positions processed')
	
	return(CHR_INFORMATIVE_POSITIONS)



def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to perform the scRNA somatic variant calling')
	parser.add_argument('--infile', type=str, help='Input file with all samples merged in a single tsv', required = True)   
	parser.add_argument('--outfile', type=str, help='Out file prefix', required = True)
	parser.add_argument('--min_cov', type=int, default = 5, help='Minimum depth of coverage to consider a sample. [Default: 5]', required = False)
	parser.add_argument('--min_cells', type=int, default = 5, help='Minimum number of cells covering a site to consider a sample. [Default: 5]', required = False)
	parser.add_argument('--min_cell_types', type=int, default = 2, help='Minimum number of cell types with enough coverage/cells at a site to be considered as a callable [Default: 2]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	min_cov = args.min_cov
	min_cells = args.min_cells
	min_cell_types = args.min_cell_types

	# How to split genome (window sizes)
	window = 20000

	# 1. Variant calling
	print ('\n------------------------------')
	print ('Variant calling')
	print ('------------------------------\n')

	# 1.1: Step 1: Beta-binomial test and other filters
	print ('- Variant calling step 1\n')
	outfile1 = outfile + ".calling.step1.tsv"
	start = timeit.default_timer()
	informative_positions = callable_sites(infile,min_cells,min_cov,min_cell_types)

	# Coverage report
	outfile_report_all = outfile + '.coverage_cell_count.report.tsv'
	outfile_report_all_chrom = outfile + '.coverage_cell_count.per_chromosome.report.tsv'
	
	# Coverages per choromosome
	count_df = 0
	for chrom in informative_positions.keys():
		informative_positions_temp = informative_positions[chrom]
		DF_temp = pd.DataFrame.from_dict({(i,j): informative_positions_temp[i][j] 
								for i in informative_positions_temp.keys() 
								for j in informative_positions_temp[i].keys()},
								orient='index').reset_index()
		DF_temp = DF_temp.rename(columns={'level_0': 'Cell_type', 'level_1' : 'Count_category'})
		DF_temp.insert(0, '#CHROM', chrom)
		if (count_df == 0):
			DF = DF_temp
			count_df = count_df + 1
		else:
			DF = DF.append(DF_temp, ignore_index=True)

	DF.to_csv(outfile_report_all_chrom, index=False, sep = '\t')

	# Total coverage per cell type
	f = {x:'sum' for x in DF.columns if not x in ['#CHROM','Cell_type','Count_category']}
	DFtotal = DF.groupby(['Cell_type','Count_category']).agg(f).reset_index()
	DFtotal.to_csv(outfile_report_all, index=False, sep = '\t')

	# # 1.2: Step 2: Add distance, editing and PoN filters
	# print ('\n- Variant calling step 2\n')
	# print("	> Editing file used: " , editing)
	# print("	> PoN file used: " , pon)
	# outfile2 = outfile + '.somatic.calling.tsv'
	# variant_calling_step2(outfile1,distance,editing,pon,dict_variants,window,outfile2)
	# os.remove(outfile1)

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')
