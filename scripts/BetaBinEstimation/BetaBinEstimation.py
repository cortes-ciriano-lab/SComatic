#!/bin/python3
from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.rinterface import RRuntimeWarning
import os
import sys
import random
import argparse
import pandas as pd
import timeit

def estimate_read_size_and_number_of_lines(filename,out_temp):
	
	filesizeBytes = os.path.getsize(filename)

	with open(filename,'r') as f, open(out_temp,'w') as f2:
		count = 0
		for line in f:
			count = count + 1

			f2.write(line)

			if count >= 10000:
				break

	filesizeBytes2 = os.path.getsize(out_temp)
	bytesPerLine = filesizeBytes2/count
	os.remove(out_temp)

	totalLinesEst = filesizeBytes / bytesPerLine
	return(bytesPerLine,totalLinesEst)

def random_selection(file, n, bytesPerLine, totalLinesEst,out):
	# exact output count
	resultlines = 0

	# Max number of tries
	# To avoid to print many times the same line or an infinite loop
	max_tries = 0

	# Move across the file
	with open(file,'r') as fh:
		while resultlines < int(n) and max_tries < 10:
			linestart = random.randint(0, int(totalLinesEst))
			readstart = (linestart * bytesPerLine) - bytesPerLine
			if readstart < bytesPerLine:
				readstart = 0
			else:
				fh.seek(readstart)
			scratch = fh.readline()
			line = fh.readline()

			#print (line)
			if not line.startswith('#') and not line == '':
				out.write(line)
				resultlines += 1
			else:
				max_tries = max_tries + 1



def extract_sites_for_betabinom_estimation(file):

	Alleles = ["A","C","T","G","I","D","N","O"]

	CCounts = []

	with open(file, 'r') as f:
		for line in f:
			if not line.startswith('#'):
				line = line.rstrip('\n')
				
				# Create dictionary with cell types and indexes for each one
				elements = line.split('\t')

				# Reference allele
				REF = elements[2]

				# Column with read/cell count info
				INFO_i = elements[4]
				
				if not INFO_i.startswith('NA'):
					DP,NC,CC,BC,BQ,BCf,BCr = INFO_i.split('|')
					
					cc = CC.split(":")
					bc = BC.split(":")

					ALT1 = [int(cc[x]) for x in range(len(cc)) if Alleles[x] != REF]
					ALT2 = [int(bc[x]) for x in range(len(bc)) if Alleles[x] != REF]


					# Getting alt read counts
					Alt_BC = sum(ALT2)
					Ref_BC = int(DP) - Alt_BC

					# Getting alt cell counts
					Alt_CC = sum(ALT1)
					Ref_CC = int(NC) - Alt_CC
					PROP_CELLS_alt = Alt_CC/float(NC)
					PROP_READS_alt = Alt_BC/float(DP)

					if PROP_CELLS_alt < 0.15 and PROP_READS_alt < 0.15:
						# It is faster if we process data by windows
						VALUES = [Alt_CC, Ref_CC, Alt_BC, Ref_BC]
						#print ('Value\t' + '\t'.join([str(x) for x in VALUES]))
						CCounts.append(VALUES)

	return (CCounts)

def betabin_estimation(CCounts):
	
	try:
		BETABIN_params = robjects.r('''
			beta_binom_params <- function (x, y){
			  suppressPackageStartupMessages(library(VGAM))
			  fit = invisible(vglm(cbind(x,y) ~ 1, betabinomialff, trace = FALSE))
			  alpha <- Coef(fit)[[1]]
			  beta <- Coef(fit)[[2]]
			  return(c(alpha,beta))
			}
			''')

		# Sample sites from each cell type
		Alt_CC= [int(x[0]) for x in CCounts]
		Ref_CC = [int(x[1]) for x in CCounts]
		Alt_BC = [int(x[2]) for x in CCounts]
		Ref_BC = [int(x[3]) for x in CCounts]	

		# Estimation for cell counts
		PARAMS = BETABIN_params(robjects.IntVector(Alt_CC), robjects.IntVector(Ref_CC))
		alpha1 = float(tuple(PARAMS)[0])
		beta1 = float(tuple(PARAMS)[1])

		# Estimation for cell counts
		PARAMS2 = BETABIN_params(robjects.IntVector(Alt_BC), robjects.IntVector(Ref_BC))
		alpha2 = float(tuple(PARAMS2)[0])
		beta2 = float(tuple(PARAMS2)[1])
		print ('\n- Beta-binomial estimation finished:')
		print ('  > Alpha 1 (CC):', str(alpha1))
		print ('  > Beta 1 (CC):', str(beta1))
		print ('  > Alpha 2 (BC):', str(alpha2))
		print ('  > Beta 2 (BC):', str(beta2))

		return(alpha1,beta1,alpha2,beta2)

	except:
		print ('\n- Beta-binomial estimation failed.')
		sys.exit()



def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to estimate the Beta-binomial distribution parameters (alpha and beta) to be used afterwards in the BaseCellCalling.step1.py')
	parser.add_argument('--in_tsv', type=str, help='File listing the tsv files to be used for the beta-binomial fitting (obtained with BaseCellCounter.py script)', required = True)   
	parser.add_argument('--outfile', type=str, help='Report with the estimated Beta-binomial parameters', required = True)
	parser.add_argument('--n_sites', type=int, default = 500000, help='Approximate number of sites to be used for fitting the Beta-binomial distribution [Default: 500000]', required = False)
	parser.add_argument('--seed', type=int, default = 1992, help='Random seed for computation [Default: 1992]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	in_tsv = args.in_tsv
	outfile = args.outfile
	n_sites = args.n_sites
	seed = args.seed


	# 1. Random seed
	random.seed(seed)

	# 2. Temp files
	out_temp1 = outfile + '.temp1'
	out_temp2 = outfile + '.temp2'

	# 3. Load files names (if exist) in a list
	tsv_files_list = []
	with open(in_tsv, 'r') as f:
		for tsv in f:
			tsv = tsv.rstrip('\n')
			if os.path.isfile(tsv):
				tsv_files_list.append(tsv)

	print ('\n- Number of tsv files to be used:')
	print ('  > ', str(len(tsv_files_list)))

	# 4. Sample equal number of sites for each tsv
	sample_n = round(n_sites / len(tsv_files_list))

	print ('\n- Random selection of sites:')
	print ('  > Total sites:', str(n_sites))
	print ('  > Sites per sample (rounded):', str(sample_n))

	with open(out_temp2,'w') as f2:
		for file in tsv_files_list:
			# Fast and random selection of lines based on estimation of the file size
			bytesPerLine,totalLinesEst = estimate_read_size_and_number_of_lines(file,out_temp1)
			random_selection(file, sample_n, bytesPerLine, totalLinesEst,f2)

	# 5. Get values for beta-bin estimation in the temporary file created in previous step
	COUNTS = extract_sites_for_betabinom_estimation(out_temp2)
	os.remove(out_temp2)

	# 6. Estimate alpha and beta parameters for beta-bin distribution
	print ('\n- Fitting beta-binomial estimation...')
	alpha1,beta1,alpha2,beta2 = betabin_estimation(COUNTS)

	# 7. Save resuts in a mini report
	beta_binomial_parameters = pd.DataFrame([[alpha1,beta1, alpha2, beta2]],
		columns=['alpha1', 'beta1', 'alpha2','beta2'])

	beta_binomial_parameters.to_csv(outfile,index=False, sep = '\t')

	print ('\n- Beta-binomial estimation report saved in:')
	print ('  > ', str(outfile))

# -------------------------
# Beta binomial
# -------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')
