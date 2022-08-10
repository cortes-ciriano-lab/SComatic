#!/usr/bin/python

import timeit
import os
import argparse
import subprocess
import time
import itertools
import pandas as pd
import numpy as np

def longestRun(s):
	if len(s) == 0: return 0
	runs = ''.join('*' if x == y else ' ' for x,y in zip(s,s[1:]))
	starStrings = runs.split()
	if len(starStrings) == 0: return 1
	return 1 + max(len(stars) for stars in starStrings)

def trincucleotide_geneartor(l):
	yield from itertools.product(*([l] * 3)) 

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to obtain trincucleotide context background')
	parser.add_argument('--in_tsv', type=str, help='File listing the tsv files to be used for the trinucleotide context background computation (files obtained in BaseCellCalling.step1.py)', required = True)   
	parser.add_argument('--out_file', type=str, help='Output file', required = True)
	return (parser)

def main():

	# 0. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	in_tsv = args.in_tsv
	outfile = args.out_file

	# in_tsv = '/nfs/research/icortes/DATA/scRNA_somatic_variant_calling/DATA/Li_Cell2020_cSCC/results/step1.targeted_regions.files.tsv'
	# outfile = '/nfs/research/icortes/DATA/scRNA_somatic_variant_calling/DATA/Li_Cell2020_cSCC/results/test.tsv'

	# Temp file and folder
	temp = outfile + '.temp'
	temp_folder=os.path.dirname(outfile)

	# 1. Concatenate and create temp calling file
	print ('-----------------------------------------------------------')
	print ('1. Preparing temp file...')
	print ('-----------------------------------------------------------\n')

	command = "for tsv in $(cat %s);do grep -v '^#\\|LC_' $tsv | awk -F'\\t' -v OFS='\\t' '{if ($19 >= 2 && $20 >= 2) {print $4,$8,$9}}'; done | sort -T %s | uniq -c | sed 's/^[ \\t]*//' | tr ' ' '\\t' > %s" % (in_tsv,temp_folder,temp)
	# Submit linux command
	try:
		subprocess.run(command, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)

	# 2. Checking trinucleotide context background 
	print ('-----------------------------------------------------------')
	print ('2. Checking trinucleotide context background...')
	print ('-----------------------------------------------------------\n')

	# Complementary bases
	complemetary_dictionary = {'A':'T','T':'A','G':'C','C':'G'}
	
	# Expected bases
	expected_bases = ['C','T']

	# Generate all possible background reference trinucleotides
	all_trinucleotides = {}
	for x in trincucleotide_geneartor('ACTG'):
		SEQ=''.join(x)
		
		if SEQ[1] in expected_bases:
			all_trinucleotides[SEQ] = 0

	# Checking and counting trinucleotides
	with open(temp,'r') as t:
		for line in t:
			line = line.rstrip('\n')
			COUNT, REF, UP, DOWN = line.split("\t")
			COUNT = int(COUNT)

			print (line)
			
			#print (UP, DOWN)
			try:
				if (longestRun(UP) < 4 and longestRun(DOWN) < 4):
					UP_base = UP[-1]
					DOWN_base = DOWN[0]

					if (REF in 'ACTG' and UP_base in 'ACTG' and DOWN_base in 'ACTG'):
						TRINUCLEOTIDE = [UP_base, REF, DOWN_base]

						if REF in expected_bases:
							trincucleotide = ''.join(TRINUCLEOTIDE)
						else:
							TRINUCLEOTIDE2 = [complemetary_dictionary[x] for x in TRINUCLEOTIDE]
							trincucleotide = ''.join(TRINUCLEOTIDE2)

						all_trinucleotides[trincucleotide] = all_trinucleotides[trincucleotide] + COUNT
			except:
				continue
	
	# Transform to data frame	
	df = pd.DataFrame.from_dict(all_trinucleotides, orient='index')
	df.index.name = 'Trincucleotide'
	df.reset_index(inplace=True)
	SUM = (np.sum(df[0]))
	df['Ratio'] = df[0] / float(SUM)
	df.rename({0: 'Count'}, axis=1, inplace=True)
	df.to_csv(outfile, index=False)

	# # Remove temp file
	# os.remove(temp)

# -------------------------
# PoN construction
# -------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')




