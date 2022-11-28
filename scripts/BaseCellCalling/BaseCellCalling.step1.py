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

def variant_calling_step1(file,alpha1,beta1,alpha2,beta2,min_ac_cells,min_ac_reads,min_cells,min_reads,min_cell_types,max_cell_types,Fisher_cutoff, window,outfile,FASTA):
	Alleles = ["A","C","T","G","I","D","N","O"]
	outfile= open(outfile,'w')
	VARIANT_POSITIONS = {}
	CHR_INFORMATIVE_POSITIONS = {}

	cur_chr = 'z'

	# Load reference file. Mandatory to be done inside function to avoid overlap problems during multiprocessing
	if FASTA != None:
		inFasta = pysam.FastaFile(FASTA)

	counter = 0
	with open(file, 'r') as f:
		for line in f:
			if line.startswith('##'):
				outfile.write(line)

			else:
				if line.startswith('#CHROM') and counter == 0:
					line = line.rstrip('\n')

					# Create dictionary with cell types and indexes for each one
					elements = line.split('\t')
					cell_types_idx = {x : elements[x] for x in range(len(elements)) if x > 4} # First 4 columns are common across all lines ['#CHROM', 'Start','End', 'REF', 'INFO']. The rest are cell types

					diff_cell_types = list(cell_types_idx.values())
					
					# Info header
					DICT = {'ALT': "##INFO=ALT,Description=Alternative alleles found",
						'FILTER':"##INFO=FILTER,Description=Filter status of the variant site",
						'Cell_types':"##INFO=Cell_types,Description=Cell type/s with the variant",
						'Up_context':"##INFO=Up_context,Description=Up-stream bases in reference (4 bases)",
						'Down_context':"##INFO=Down_context,Description=Down-stream bases in reference (4 bases)",
						'N_ALT':"##INFO=N_ALT,Description=Cell type/s with the variant",
						'Dp':"##INFO=Dp,Description=Depth of coverage (reads) in the cell type supporting the variant",
						'Nc':"##INFO=Nc,Description=Number of distinct cells found in the cell type with the mutation",
						'Bc':"##INFO=Bc,Description=Number of reads (base count) supporting the variants in the cell type with the mutation",
						'Cc':"##INFO=Cc,Description=Number of distinct cells supporting the variant in the cell type with the mutation",
						'VAF':"##INFO=VAF,Description=Variant allele frequency of variant in the cell type with the mutation",
						'CCF':"##INFO=CCF,Description=Cancer cell fraction (fraction of ditinct cells) supporting the alternative allele in the cell type with the mutation",
						'BCp':"##INFO=BCp,Description=Beta-binomial p-value for the variant allele (considering read counts)",
						'CCp':"##INFO=CCp,Description=Beta-binomial p-value for the variant allele (considering cell counts)",
						'Cell_types_min_BC':"##INFO=Cell_types_min_BC,Description=Number of cell types with a minimum number of reads covering a site",
						'Cell_types_min_CC':"##INFO=Cell_types_min_CC,Description=Number of cell types with a minimum number of distinct cells found in a specific site",
						'Rest_BC':"##INFO=Rest_BC,Description=Base counts (reads) supporting other alternative alleles in this site. BC;DP;P-value (betabin)",
						'Rest_CC':"##INFO=Rest_CC,Description=Cell counts supporting other alternative alleles in this site. CC;NC;P-value (betabin)",
						'Fisher_p':"##INFO=Fisher_p,Description=Strand bias test. Fisher exact test p-value between forward and reverse reads in variant and reference allele",
						'Cell_type_Filter':"##INFO=Cell_type_Filter,Description=Filter status of the variant site in each cell type"
						}
					for info in DICT.keys():
						outfile.write(DICT[info]+'\n')	

					# Add calling columns
					INFO = ['ALT', 'FILTER', 'Cell_types','Up_context', 'Down_context','N_ALT', 'Dp', 'Nc', 'Bc', 'Cc', 'VAF', 'CCF', 'BCp', 'CCp', 'Cell_types_min_BC', 'Cell_types_min_CC', "Rest_BC", "Rest_CC",'Fisher_p', 'Cell_type_Filter']
					CALLING_INFO = "\t".join(INFO)
					elements.insert(4,CALLING_INFO)
					outfile.write('\t'.join(elements)+'\n')

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

					# Context regions (5 bases up- and down-stream)
					num_bases = 5
					if FASTA != None:
						try:
							context = inFasta.fetch(CHROM, int(POS)-(num_bases+1), int(POS)+num_bases)
							context = context.upper()
							up_context = context[0:num_bases]
							down_context = context[(num_bases+1):(num_bases+1+num_bases)]
						except:
							up_context = '.'
							down_context = '.'
					else:
						up_context = '.'
						down_context = '.'


					# Reference allele
					REF = elements[3]

					### Single cell type level info
					## Variant calling info necessary to make a call
					# Significant alternative alleles
					Alts = []
					# Significant cell types
					Cell_types = []
					# Coverage of the significant cell types
					DPs = []
					# Number of cells in significant cell types
					NCs = []
					# Base counts of the significant alleles
					BCs = []
					# Cell counts of the significant alleles
					CCs = []
					# Beta-bin p-value of the significant alleles (base counts)
					BCp = []
					# Beta-bin p-value of the significant alleles (cell counts)
					CCp = []
					# Variant allele frequency (base on read counts)
					VAF = []
					# Cancer cell fraction (VAF base on cell counts)
					CCF = []
					# Number of cell types with min_BC
					Cell_types_min_BC = 0
					# Number of cell types with min_CC
					Cell_types_min_CC = 0
					# Filter
					Filter = []
					# Fisher strand
					Fisher_p = []
					# Sum alts
					Sum_alts_bc = 0
					# Sum alts cc
					Sum_alts_cc = 0
					# Dp_sums 
					Sum_dp = 0
					# Sum number of cells
					Sum_nc = 0
					# All DPs
					All_DPs = {}
					# All CCs
					All_CCs = {}

					### Site variant calling level info
					FILTER = []

					# Computation for each cell type
					for cell_type_i in cell_types_idx.keys():
						cell_type = cell_types_idx[cell_type_i]
						INFO_i = elements[cell_type_i]

						if INFO_i.startswith('NA'):
							All_DPs[cell_type] = 0
							All_CCs[cell_type] = 0
 
						elif not INFO_i.startswith('NA'):
							DP,NC,CC,BC,BQ,BCf,BCr = INFO_i.split('|')
							
							All_DPs[cell_type] = int(DP)
							All_CCs[cell_type] = int(NC)

							if (int(DP) >= min_cells and int(NC) >= min_reads):
								# Number of cell types with min_BC
								Cell_types_min_BC = Cell_types_min_BC + 1
								# Number of cell types with min_CC
								Cell_types_min_CC = Cell_types_min_CC + 1
								
								cc = CC.split(":")
								bc = BC.split(":")
								bq = BQ.split(":")
								bcf = BCf.split(":") # forward bases
								bcr = BCr.split(":") # reverese bases

								# Sum all alternative reads
								Alts2 = sum([int(bc[x]) for x in range(len(bc)) if Alleles[x] not in [REF,"O"]])
								Sum_alts_bc = Sum_alts_bc + Alts2
								Cc2 = sum([int(cc[x]) for x in range(len(cc)) if Alleles[x] not in [REF,"O"]])
								Sum_alts_cc = Sum_alts_cc +Cc2
								Sum_dp = Sum_dp + int(DP)
								Sum_nc = Sum_nc + int(NC)

								# Get p-values for base counts
								Alt_bc_dict = {Alleles[x]:int(bc[x]) for x in range(len(bc)) if Alleles[x] not in [REF,"I","D","N","O"] and int(bc[x]) > 0}
								Alt_bc_p_dict = {x : round(betabinom.sf(Alt_bc_dict[x]-0.1, int(DP), alpha1, beta1),4) for x in Alt_bc_dict.keys() }
								Alt_bc_p_dict_sig = [x for x in Alt_bc_p_dict.keys() if Alt_bc_p_dict[x] < 0.001]

								# Get p-values for cell counts
								Alt_cc_dict = {Alleles[x]:int(cc[x]) for x in range(len(cc)) if Alleles[x] not in [REF,"I","D","N","O"] and int(cc[x]) > 0}
								Alt_cc_p_dict = {x : round(betabinom.sf(Alt_cc_dict[x]-0.1, int(NC), alpha2, beta2),4) for x in Alt_cc_dict.keys() }
								Alt_cc_p_dict_sig = [x for x in Alt_cc_p_dict.keys() if Alt_cc_p_dict[x] < 0.001]

								Alt_candidates = sorted(list(set(Alt_bc_p_dict_sig + Alt_cc_p_dict_sig)))

								if (len(Alt_candidates) > 0):
									# Updating required lists for calling
									alt_candidates = "|".join(Alt_candidates)
									Alts.append(alt_candidates)
									
									Cell_types.append(cell_type)
									
									DPs.append(str(DP))
									
									NCs.append(str(NC))

									# Probabilities of Base counts
									P_BC = list(Alt_bc_p_dict.values())

									# Probabilities of Cell counts
									P_CC = list(Alt_cc_p_dict.values())

									# Strand bias test (Fisher's exact test between strands, comparing alt and ref distributions)
									if (Fisher_cutoff != 1):
										Fw_dict = {Alleles[x]:int(bcf[x]) for x in range(len(bcf)) if Alleles[x]}
										Rv_dict = {Alleles[x]:int(bcf[x]) for x in range(len(bcf)) if Alleles[x]}
										fisher_p = "|".join([str(round(stats.fisher_exact([[Fw_dict[x], Rv_dict[x]], [Fw_dict[REF], Rv_dict[REF]]])[1],4)) for x in Alt_candidates])
										Fisher_p.append(fisher_p)

									b = "|".join([str(Alt_bc_dict[x]) for x in Alt_candidates])
									BCs.append(b)
									
									c = "|".join([str(Alt_cc_dict[x]) for x in Alt_candidates])
									CCs.append(c)

									bp = "|".join([str(Alt_bc_p_dict[x]) for x in Alt_candidates])
									BCp.append(bp)

									cp = "|".join([str(Alt_cc_p_dict[x]) for x in Alt_candidates])
									CCp.append(cp)

									vaf = "|".join([str(round(Alt_bc_dict[x]/float(DP),4)) for x in Alt_candidates])
									VAF.append(vaf)

									ccf = "|".join([str(round(Alt_cc_dict[x]/float(NC),4)) for x in Alt_candidates])
									CCF.append(ccf)

									# Save information about not candidate cells
									b0 = sum([int(Alt_bc_dict[x]) for x in Alt_candidates])
									c0 = sum([int(Alt_cc_dict[x]) for x in Alt_candidates])
									Sum_dp = Sum_dp - b0
									Sum_nc = Sum_nc - c0
									Sum_alts_bc = Sum_alts_bc - b0
									Sum_alts_cc = Sum_alts_cc - c0

									# Check if position passes or not
									# Cell type variant calling
									# 1. Check if both read counts (base counts) and cell counts passed the p-value filters
									if not (min(P_BC) < 0.001 and min(P_CC) < 0.001):
										Filter.append('BetaBin_problem')
									elif len(Alt_candidates) > 1: # Check for multiallelic sites
										Filter.append('Multi-allelic')
									elif int(c) < min_ac_cells: # Check number of cells supporting alternative allele
										Filter.append('Low_cells')
									elif int(b) < min_ac_reads: # Check number of reads supporting alternative allele
										Filter.append('Low_reads')
									elif Fisher_cutoff != 1: # Check fisher p-value
										if float(fisher_p) < Fisher_cutoff:
											Filter.append('Fisher')
									else:
										Filter.append('PASS')


					# Save coverage and cell count information
					if (Cell_types_min_CC >= min_cell_types):
						# DPs
						for cell_type in All_DPs.keys():
							DP = All_DPs[cell_type]
							NC = All_CCs[cell_type]

							for x in TEMPLATE.keys():
								if DP >= x:
									CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['DP'][x] = CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['DP'][x] + 1
								if NC >= x: 
									CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['NC'][x] = CHR_INFORMATIVE_POSITIONS[CHROM][cell_type]['NC'][x] + 1

					# 2. Site variant calling (considering all cell types into account)
					if (len(Alts) > 0):
						# Store variant positions to check the distance between variants later on
						WIND = math.floor(POS / float(window))
						if CHROM in VARIANT_POSITIONS.keys():
							if WIND in VARIANT_POSITIONS[CHROM].keys():
								VARIANT_POSITIONS[CHROM][WIND].append(POS)
							else:
								VARIANT_POSITIONS[CHROM][WIND] = [POS]
						else:
							VARIANT_POSITIONS[CHROM] = {}
							VARIANT_POSITIONS[CHROM][WIND] = [POS]

						# 1. Filter for how many pass cell types we have
						PASS = [x for x in Filter if x == 'PASS']
						if (len(PASS) > max_cell_types):
							FILTER.append('Multiple_cell_types')

						# 2. Diffent significant alleles found
						LEN_Alts = len(set(Alts))
						if (LEN_Alts > 1 or 'Multi-allelic' in Filter):
							FILTER.append('Multi-allelic')

						# 3. Minimum number of cell types with enough info
						if (Cell_types_min_CC < min_cell_types):
							FILTER.append('Min_cell_types')

						# 4. Cell type noise
						if (len(Filter) - len(PASS) > 0):
							FILTER.append('Cell_type_noise')

						# 5. Noisy site
						# This filter check if, when ignore cells or reads harbouring the candidate mutations, there is significant noise in the site
						# It collapses all reads and cells across different cell types
						if (Sum_alts_bc > 0):
							BC_noise_p = round(1 - betabinom.cdf(Sum_alts_bc-0.1, Sum_dp, alpha1, beta1),4)
							CC_noise_p = round(1 - betabinom.cdf(Sum_alts_cc-0.1, Sum_nc, alpha2, beta2),4)
							rest_BC = [str(Sum_alts_bc), str(Sum_dp),str(BC_noise_p)]
							rest_CC = [str(Sum_alts_cc), str(Sum_nc),str(CC_noise_p)]
						else:
							BC_noise_p = 1
							CC_noise_p = 1
							rest_BC = [str(Sum_alts_bc), str(Sum_dp),str(BC_noise_p)]
							rest_CC = [str(Sum_alts_cc), str(Sum_nc),str(CC_noise_p)]

						rest_BC = ";".join(rest_BC)
						rest_CC = ";".join(rest_CC)

						if (BC_noise_p < 0.05 or CC_noise_p < 0.05):
							FILTER.append('Noisy_site')

						# 6. Context filter
						# Up-stream
						homopolymer_up = homopolymer_function(up_context,Alts,'upstream')
						if homopolymer_up == 1:
							FILTER.append("LC_Upstream")

						# Down-stream
						homopolymer_down = homopolymer_function(down_context,Alts,'downstream')
						if homopolymer_down == 1:
							FILTER.append("LC_Downstream")

						# Final filter. Collapse all FILTER information
						if (len(FILTER) == 0):
							if 'PASS' in Filter:
								FILTER = 'PASS'
							else:
								FILTER = ",".join(Filter)
						else:
							FILTER = ",".join(FILTER)

						# 6. Get calling columns for each candidate site
						Alts = ",".join(Alts)
						# Significant cell types
						Cell_types = ",".join(Cell_types)
						# Coverage of the significant cell types
						DPs = ",".join(DPs)
						# Number of cells in significant cell types
						NCs = ",".join(NCs)
						# Base counts of the significant alleles
						BCs = ",".join(BCs)
						# Cell counts of the significant alleles
						CCs = ",".join(CCs)
						# VAF
						VAF = ",".join(VAF)
						# CCF
						CCF = ",".join(CCF)
						# Beta-bin p-value of the significant alleles (base counts)
						BCp = ",".join(BCp)
						# Beta-bin p-value of the significant alleles (cell counts)
						CCp = ",".join(CCp)
						# Fisher
						if (Fisher_cutoff != 1):
							Fisher_p = ",".join(Fisher_p)
						else:
							Fisher_p = '.'
						# Number of cell types with min_BC
						Cell_types_min_BC = str(Cell_types_min_BC)
						# Number of cell types with min_CC
						Cell_types_min_CC = str(Cell_types_min_CC)
						# Filter
						Filter = ",".join(Filter)

						# Create out line with call performed
						INFO = [Alts, FILTER, Cell_types, up_context, down_context, str(LEN_Alts), DPs, NCs, BCs, CCs, VAF, CCF, BCp, CCp, Cell_types_min_BC, Cell_types_min_CC, rest_BC, rest_CC, Fisher_p, Filter]
						CALLING_INFO = "\t".join(INFO)
						elements.insert(4,CALLING_INFO)
						outfile.write('\t'.join(elements)+'\n')
					else:
						Alts = "."
						LEN_Alts = '.'
						# Significant cell types
						Cell_types = "."
						# Coverage of the significant cell types
						DPs = "."
						# Number of cells in significant cell types
						NCs = "."
						# Base counts of the significant alleles
						BCs = "."
						# Cell counts of the significant alleles
						CCs = "."
						# Beta-bin p-value of the significant alleles (base counts)
						BCp = "."
						# Beta-bin p-value of the significant alleles (cell counts)
						CCp = "."
						# Variant allele frequency (base on read counts)
						VAF = "."
						# Cancer cell fraction (VAF base on cell counts)
						CCF = "."
						# Fisher
						Fisher_p = '.'
						# Rest filter
						if (Sum_alts_bc > 0):
							BC_noise_p = round(1 - betabinom.cdf(Sum_alts_bc-0.1, Sum_dp, alpha1, beta1),4)
							CC_noise_p = round(1 - betabinom.cdf(Sum_alts_cc-0.1, Sum_nc, alpha2, beta2),4)
							rest_BC = [str(Sum_alts_bc), str(Sum_dp),str(BC_noise_p)]
							rest_CC = [str(Sum_alts_cc), str(Sum_nc),str(CC_noise_p)]
						else:
							BC_noise_p = 1
							CC_noise_p = 1
							rest_BC = [str(Sum_alts_bc), str(Sum_dp),str(BC_noise_p)]
							rest_CC = [str(Sum_alts_cc), str(Sum_nc),str(CC_noise_p)]

						rest_BC = ";".join(rest_BC)
						rest_CC = ";".join(rest_CC)

						if (BC_noise_p < 0.001 or CC_noise_p < 0.001):
							Filter = '.'
							FILTER = 'Noisy_site' 

							# Save position of noise in the list
							WIND = math.floor(POS / float(window))
							if CHROM in VARIANT_POSITIONS.keys():
								if WIND in VARIANT_POSITIONS[CHROM].keys():
									VARIANT_POSITIONS[CHROM][WIND].append(POS)
								else:
									VARIANT_POSITIONS[CHROM][WIND] = [POS]
							else:
								VARIANT_POSITIONS[CHROM] = {}
								VARIANT_POSITIONS[CHROM][WIND] = [POS]
						else:
							Filter = "."
							FILTER = '.'

						# Number of cell types with min_BC
						Cell_types_min_BC = str(Cell_types_min_BC)
						# Number of cell types with min_CC
						Cell_types_min_CC = str(Cell_types_min_CC)

						# Create out line with call performed
						INFO = [Alts, FILTER, Cell_types, up_context, down_context, str(LEN_Alts), DPs, NCs, BCs, CCs, VAF, CCF, BCp, CCp, Cell_types_min_BC, Cell_types_min_BC, rest_BC, rest_CC,Fisher_p, Filter]
						CALLING_INFO = "\t".join(INFO)
						elements.insert(4,CALLING_INFO)
						outfile.write('\t'.join(elements)+'\n')

					if (counter % 1000000 == 0):
						print ('Step 1 variant calling: ' + str(counter) + ' positions processed')
	
	if FASTA != None:
		inFasta.close()

	outfile.close()
	return(VARIANT_POSITIONS,CHR_INFORMATIVE_POSITIONS)

def longestRun(s):
	if len(s) == 0: return 0
	runs = ''.join('*' if x == y else ' ' for x,y in zip(s,s[1:]))
	starStrings = runs.split()
	if len(starStrings) == 0: return 1
	return 1 + max(len(stars) for stars in starStrings)

def FrequentBase(s):
	L = list(s)
	
	# Get the most common element and the percentage respect the length of the sequence
	MAX = max([L.count(base) for base in set(L)])
	PERC = round(float(MAX)/len(s), 2)
	
	return(PERC)

def homopolymer_function(A):
	# A is the up or down-stream sequence
	if (A != '.'):
		# Get the longest k-mer
		Max_seq_l = longestRun(A)
		# Get the frequency of the more repeated base
		Max_base_freq = FrequentBase(A)
		
		if (Max_seq_l > 3 and Max_base_freq >= 0.8):
			FILTER = 1
		else:
			FILTER = 0
	else:
		FILTER = 0
		
	return(FILTER)

def homopolymer_function(A,alts,stream_direction):
	# A is the up or down-stream sequence
	if (A != '.'):

		# Get the longest k-mer including the alt base
		if (stream_direction == 'upstream'):
			Max_seq_l = max([longestRun(A+x) for x in alts])
		elif (stream_direction == 'downstream'):
			Max_seq_l = max([longestRun(x+A) for x in alts])
		
		# If the longest homo-kmer (including the alternative allele) is longer than 4, filter this site
		if (Max_seq_l >= 4):
			FILTER = 1
		else:
			FILTER = 0
	else:
		FILTER = 0
		
	return(FILTER)

def GetClosestHetSites(CHROM, SITE, HET_DICT, window,distance):
	# Empty variables (just in case)
	LIST_HET = []
	LIST_HET1 = []

	# Get chromosome window for candidate site +/- the nightbour one 
	TARGET_WIND_SITE = math.floor(SITE / float(window))
	WIND_SITE_LIST = [TARGET_WIND_SITE - 1, TARGET_WIND_SITE, TARGET_WIND_SITE + 1]

	# Find the Het sites
	LIST_HET = []
	for WIND_SITE in WIND_SITE_LIST:
		try:
			LIST_HET0 = HET_DICT[CHROM][WIND_SITE]
			LIST_HET1 = [x for x in LIST_HET0 if (abs(int(x) - SITE) <= distance and abs(int(x) - SITE) > 0)]
		except:
			LIST_HET1 = []

		# Concatenate new lists
		LIST_HET = LIST_HET + LIST_HET1
		if len(LIST_HET) > 0:
			return(1)
	
	# If het list is too long, we only select the closest 100 het sites
	return(len(LIST_HET))


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
	parser.add_argument('--outfile', type=str, help='Output file prefix', required = True)
	parser.add_argument('--ref', type=str, help='Reference fasta file (*fai must exist)', required = True)
	parser.add_argument('--editing', type=str, help='RNA editing file to be used to remove RNA-diting sites', required = False)
	parser.add_argument('--pon', type=str, help='Panel of normals (PoN) file to be used to remove germline and false positive calls', required = False)	
	parser.add_argument('--min_cov', type=int, default = 5, help='Minimum depth of coverage to consider a sample. [Default: 5]', required = False)
	parser.add_argument('--min_cells', type=int, default = 5, help='Minimum number of cells covering a site to consider a sample. [Default: 5]', required = False)
	parser.add_argument('--min_ac_cells', type=int, default = 2, help='Minimum number of cells supporting the alternative allele to consider a mutation. [Default: 2]', required = False)
	parser.add_argument('--min_ac_reads', type=int, default = 3, help='Minimum number of reads supporting the alternative allele to consider a mutation. [Default: 3]', required = False)
	parser.add_argument('--max_cell_types', type=int, default = 1, help='Maximum number of cell types carrying a mutation to make a somatic call. [Default: 1]', required = False)
	parser.add_argument('--min_cell_types', type=int, default = 2, help='Minimum number of cell types with enough coverage and cell to consider a site as callable [Default: 2]', required = False)
	parser.add_argument('--fisher_cutoff', type=float, default = 1, help='P-value cutoff for the Fisher exact test performed to detect strand bias. Expected float value, if applied, we recommend 0.001. By default, this test is switched off with a value of 1 [Default: 1]', required = False)
	parser.add_argument('--min_distance', type=int, default = 5, help='Minimum distance allowed between potential somatic variants [Default: 5]', required = False)
	parser.add_argument('--alpha1', type=float, default = 0.260288007167716, help='Alpha parameter for Beta-binomial distribution of read counts. [Default: 0.260288007167716]', required = False)
	parser.add_argument('--beta1', type=float, default = 173.94711910763732, help='Beta parameter for Beta-binomial distribution of read counts. [Default: 173.94711910763732]', required = False)
	parser.add_argument('--alpha2', type=float, default = 0.08354121346569514, help='Alpha parameter for Beta-binomial distribution of cell counts. [Default: 0.08354121346569514]', required = False)
	parser.add_argument('--beta2', type=float, default = 103.47683488327257, help='Beta parameter for Beta-binomial distribution of cell counts. [Default: 103.47683488327257]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	fasta = args.ref
	editing = args.editing
	pon = args.pon
	min_cov = args.min_cov
	min_cells = args.min_cells
	min_ac_cells = args.min_ac_cells
	min_ac_reads = args.min_ac_reads
	max_cell_types = args.max_cell_types
	min_cell_types = args.min_cell_types
	fisher_cutoff = args.fisher_cutoff
	distance = args.min_distance
	alpha1 = args.alpha1
	beta1 = args.beta1
	alpha2 = args.alpha2
	beta2 = args.beta2

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
	dict_variants, informative_positions = variant_calling_step1(infile,alpha1,beta1,alpha2,beta2,min_ac_cells,min_ac_reads,min_cells,min_cov,min_cell_types,max_cell_types,fisher_cutoff,window,outfile1,fasta)

	# Coverage report
	outfile_report_all = outfile + '.coverage_cell_count.report.tsv'
	outfile_report_all_chrom = outfile + '.coverage_cell_count.per_chromosome.report.tsv'
	

#-------------------------
# Running scRNA somatic variant calling
#-------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')
