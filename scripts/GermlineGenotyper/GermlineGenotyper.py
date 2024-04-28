import pandas as pd
import pybedtools
import subprocess
import math
import argparse
import timeit
import os

#----------
# Script to extract germline information from SComatic output
#----------

'''
List of functions
'''
# Function to build a dictionary with the sites to check (SNPs)
def build_dict(snp, window):
	DICT_germ = {}
	with open(snp, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				# Coordinates and variant info
				line = line.rstrip('\n')
				ID=line
				CHROM,POS,REF,ALT = line.split(':')
				POS = int(POS)

				WIND = math.floor(POS / float(window))

				if not CHROM in DICT_germ.keys():
					DICT_germ[CHROM] = {}
					DICT_germ[CHROM][WIND] = {}
					DICT_germ[CHROM][WIND][POS] = [REF,ALT]
				else:
					if not WIND in DICT_germ[CHROM].keys():
						DICT_germ[CHROM][WIND] = {}
						DICT_germ[CHROM][WIND][POS] = [REF,ALT]
					else:
						DICT_germ[CHROM][WIND][POS] = [REF,ALT]
	return (DICT_germ)


# Function to perform the germline variant calling
def GermlineCalling(out2,out3,min_reads,min_cells,min_cov,min_total_cells,min_cell_types_exp, DICT_germ,window):
	outfile3= open(out3,'w')

	# Write header
	HEADER = ['MUT_ID','CHROM', 'POS', 'REF', 'ALT', 'GT', 'Germline_filter', 'Cell_types', 'Cell_types_min_BC', 'Cell_types_min_CC', 'Read_depth', 'Read_count_ALT', 'VAF', 'Cell_depth', 'Cell_count_ALT', 'CCF', 'SComatic_filter']
	HEADER = "\t".join(HEADER)
	outfile3.write(HEADER + '\n')

	# Run for all variants
	with open(out2, 'r') as f:
		for line in f:
			if line.startswith('#'):
				pass
				#outfile.write(line)
			else:
				line = line.rstrip('\n')
				elements = line.split('\t')

				# Important variables to consider
				CHROM = elements[0]
				POS = int(elements[2])
				REF = elements[3]
				ALT = elements[4]
				FILTER = elements[5]
				Cell_types = elements[6]
				N_ALT = elements[9]
				Dp = elements[10]
				Nc = elements[11]
				Bc = elements[12]
				Bc2 = Bc.replace('|',',')
				Cc = elements[13]
				Cc2 = Cc.replace('|',',')
				VAF = elements[14]
				CCF = elements[15]
				Cell_types_min_BC = elements[18]
				Cell_types_min_CC = elements[19]
				Rest_BC = elements[20]
				Rest_CC = elements[21]

				# Sum rest of bases
				Rest_BCs = [float(x) for x in Rest_BC.split(';')][:-1]
				Rest_ACs = int(Rest_BCs[0])
				Rest_BCs = sum(Rest_BCs)
				Rest_CCs = [float(x) for x in Rest_CC.split(';')][:-1]
				Rest_ACCs = int(Rest_CCs[0])
				Rest_CCs = sum(Rest_CCs)

				# Unique alternative alleles
				ALTs = set(ALT.split(','))
				ALTs = ",".join(ALTs)

				# Compute window and extract the expected 
				WIND = math.floor(POS / float(window))
				try:
					REF_snp,ALT_snp = DICT_germ[CHROM][WIND][POS]
				except:
					continue
				MUT_ID = f"{CHROM}:{POS}:{REF_snp}:{ALT_snp}"

				# Check if there is existing alternative allele
				# Filter column
				New_filters = []
				if (ALT == "."):
					Dps = Rest_BCs
					Ncs = Rest_CCs
					Bcs = Rest_ACs
					Ccs = Rest_ACCs
					VAFall = round(Bcs/Dps,2)
					CCFall = round(Ccs/Ncs,2)
					GT = "0/0"

					# Common filters
					if Dps < min_cov:
						New_filters.append('Low_cov')
					# Total coverage (cells)
					if Ncs < min_total_cells:
						New_filters.append('Low_total_cells')
					# Check the number of cell types carrying it
					if int(Cell_types_min_BC) < min_cell_types_exp and int(Cell_types_min_CC) < min_cell_types_exp:
						New_filters.append('Low_cell_type_exp')
				else:
					# Sum of coverage for the cell types supporting the variant
					Dps = sum([int(x) for x in Dp.split(',')])
					Dps = Dps + Rest_BCs
					# Sum of cells coveraged for the cell types supporting the variant
					Ncs = sum([int(x) for x in Nc.split(',')])
					Ncs = Ncs + Rest_CCs
					# Sum of reads with the mutations within the cell types with the variant
					Bcs = sum([int(x) for x in Bc2.split(',')])
					# Sum of cells with the mutations within the cell types with the variant
					Ccs = sum([int(x) for x in Cc2.split(',')])
					# N cell types
					N_cell_types_mutated = len(ALT.split(','))
					# All vafs
					VAFall = round(float(Bcs) / Dps,2)
					# All CCF
					CCFall = round(float(Ccs) / Ncs,2)


					# Check few filters
					# Total coverage (reads)
					if Dps < min_cov:
						New_filters.append('Low_cov')
					# Total coverage (cells)
					if Ncs < min_total_cells:
						New_filters.append('Low_total_cells')
					# Check the number of cell types carrying it
					if int(Cell_types_min_BC) < min_cell_types_exp and int(Cell_types_min_CC) < min_cell_types_exp:
						New_filters.append('Low_cell_type_exp')
					# Min Bc and Cc >= 3
					if Bcs < min_reads:
						New_filters.append('Low_read_support')
					if Ccs < min_cells:
						New_filters.append('Low_cell_support')
					# Total VAF >= 0.1
					if (VAFall < 0.1 and CCFall < 0.1):
						New_filters.append('Low_VAF_CCF')
					# Check if ALT is the expected one
					if ALT_snp not in ALTs:
						New_filters.append('ALT_not_expected')
					if int(N_ALT) > 1 or "|" in ALT:
						New_filters.append('Multi_allelic')                      

					# Genotype
					if (VAFall >= 0.1 and VAFall < 0.9):
						GT = "0/1"
					elif (VAFall >= 0.9):
						GT = "1/1"
					else:
						GT = "0/0"

				# Merge filters in a single string
				New_filters = ",".join(New_filters)
				if (New_filters == '' and GT != "0/0"):
					New_filters = "PASS"
				elif (New_filters == '' and GT == "0/0"):
					New_filters = "."

				# Final column
				INFO = [MUT_ID,CHROM, str(POS), REF, ALTs, GT, New_filters, Cell_types, str(Cell_types_min_BC), str(Cell_types_min_CC), str(Dps), str(Bcs), str(VAFall), str(Ncs), str(Ccs), str(CCFall), FILTER]
				INFO = "\t".join(INFO)
				outfile3.write(INFO + '\n')

	outfile3.close()


def GenotypeEthnicity(GT,Germline_filter):
	dict_gt = {"0/0": 0, "0/1": 0.5, "1/1": 1}
	if Germline_filter == "PASS" or Germline_filter == ".":
		GT_new = dict_gt[GT]
	else:
		GT_new = "NA"

	return GT_new

def EthnicityFormat(out3,snp,out4):

	# Mandatory rows
	SNP = pd.read_csv(snp, header = None)
	SNP.columns = ['MUT_ID']

	# Open germline calls file 
	OUT3 = pd.read_csv(out3, sep = '\t')
	OUT3 = OUT3[['MUT_ID','GT','Germline_filter']]

	# Merge results
	result = pd.merge(SNP, OUT3, on="MUT_ID", how = "left", sort = False)

	# Encode genotype for the ethnicity computation
	result['GT_new'] = result.apply(lambda row : GenotypeEthnicity(row['GT'],
						 row['Germline_filter']), axis = 1)

	# Get the column of interest
	result = result[['GT_new']]

	# Save the results
	result.to_csv(out4, sep = '\t', index = False, header = False)

'''
Arguments required for the computation of TelFusDetector
'''

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to call germline variants in a list of polymorphic sites')
	parser.add_argument('--snp', type=str, help='File with SNPs to check. Format: chr:pos:ref:alt', required = True)
	parser.add_argument('--tsv', type=str, help='TSV file obtained by the BaseCellCalling/BaseCellCalling.step1.py script', required = True)
	parser.add_argument('--prefix', default = '.', help='Prefix for the out files. Full path and prefix recommended', required = True)
	parser.add_argument('--min_reads', type=int, default = 3, help='Minimum number of reads supporting the alternative allele [Default = 3]', required = False)
	parser.add_argument('--min_cells', type=int, default = 3, help='Minimum number of cells harbouring the alternative allele [Default = 3]', required = False)
	parser.add_argument('--min_cells_types_exp', type=int, default = 1, help='Minimum number of cell types with enough expression to get a calls [Default = 1]', required = False)
	parser.add_argument('--min_cov', type=int, default = 5, help='Minimum number of reads covering the variant site. [Default = 5 ]', required = False)
	parser.add_argument('--min_total_cells', type=int, default = 5, help='Minimum number of cells with at least one read covering the variant site. [Default = 5]', required = False)
	parser.add_argument('--ethnicity', choices = ["Yes","No"], default = False, help='Specify if ethnicity format file should be created. [Default: No]', required = False)
	return (parser)


'''
Main computation
'''
def main():

	#------------
	# Get arguments
	#------------

	parser = initialize_parser()
	args = parser.parse_args()

	snp = args.snp	
	tsv = args.tsv
	prefix = args.prefix
	min_reads = args.min_reads
	min_cells = args.min_cells
	min_cov = args.min_cov
	min_total_cells = args.min_total_cells
	min_cell_types_exp = args.min_cells_types_exp


	#------------
	# 1. Preparing and processing input data
	#------------

	print ("1. Preparing and processing input data\n")

	# 0. Transform SNP file to BED and dictionary file
	out1 = f"{prefix}.snp.bed"

	SNP = pd.read_csv(snp, sep = ':', header = None)
	SNP.columns = ['CHROM','POS','REF','ALT']
	SNP['START'] = SNP['POS'] #- 1

	BED = SNP[["CHROM", "START","POS"]]
	BED.to_csv(out1, sep = '\t', index = False, header = False)

	# 1. Create a dictionary with the sites to check
	window = 50000
	DICT_germ = build_dict(snp,window)

	# 2. Intersect SComatic file with desired SNPs
	out2 = f"{prefix}.intersected.SComatic.tsv"

	# Command to run
	command = f"bedtools intersect -header -u -a {tsv} -b {out1} > {out2}" 

	# Submit linux command
	try:
		subprocess.run(command, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)

	#------------
	# 2. Running germline variant calling
	#------------

	print ("2. Running germline variant calling\n")

	# Output file
	out3 = f"{prefix}.GermlineCalls.SComatic.tsv"
	GermlineCalling(out2,out3,min_reads,min_cells,min_cov,min_total_cells,min_cell_types_exp,DICT_germ,window)

	#------------
	# 3. Prepare the output for ethnicity
	#------------

	if (args.ethnicity == "Yes"):
		print ("3. Preparing calls for the ethnicity prediction tool\n")

		# Output file
		out4 = f"{prefix}.GermlineCalls.SComatic.ethnicity.tsv"
		EthnicityFormat(out3,snp,out4)

	#------------
	# Removing temp files
	#------------
	os.remove(out1)
	os.remove(out2)

#---------------
# Running the script
#---------------

if __name__ == '__main__':
    start = timeit.default_timer()
    main()
    stop = timeit.default_timer()
    Seconds = round(stop - start)
    print(f"Computation time: {Seconds} seconds\n") 
