#!/usr/bin/python

import timeit
import os
import argparse
import subprocess
import time
# datamash in environment

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to build a SComatic Panel Of Normals (PoNs)')
	parser.add_argument('--in_tsv', type=str, help='File with tsv files to be used for final PoN construction (ideally files obtained in BaseCellCalling.step1.py)', required = True)   
	parser.add_argument('--out_file', type=str, help='PoN output file name', required = True)
	parser.add_argument('--min_samples', type=int, default = 2, help='Minimum number of significant samples to consider a site in the PoN. [Default: 2]', required = False)
	parser.add_argument('--rm_prefix', type=str, choices = ['Yes','No'], default = 'Yes', help='Remove chr prefix from input files (Yes) or no (No) [Default: Yes]', required = False)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	in_tsv = args.in_tsv
	outfile = args.out_file
	min_samples = args.min_samples
	rm_prefix = args.rm_prefix

	# 1. Merge cell types files and select potential het sites for phasing
	print ('-----------------------------------------------------------')
	print ('1. Building Panel Of Normals ...')
	print ('-----------------------------------------------------------\n')


	# Print header in outfile
	out = open(outfile,'w')
	Header=['#CHROM','POS', 'Num_samples', 'Sample_ids']
	date = time.strftime("%d/%m/%Y")## dd/mm/yyyy format
	DATE="##fileDate=%s" % date
	CONCEPTS="##INFO=Num_samples,Description=Number of significant samples (beta-binomial test)\n##INFO=Sample_ids,Description=ID of the significant samples (beta-binomial test)"

	out.write(DATE + '\n')
	out.write(CONCEPTS + '\n')
	out.write('\t'.join(Header)+'\n')
	out.close()


	# Extract recurrent sites to construct PoNs
	# Based on awk, grep, sort, and uniq unix tools

	working_dir=os.path.dirname(outfile)

	# Decide if remove prefix or not
	if (rm_prefix == 'No'): 
		command = "for file in $(cat %s);do BASE=$(basename $file); grep -v \"^#\" $file | awk -F'\\t' -v OFS='\\t' -v var=\"$BASE\" '{if ($6 != \".\") {print $1,$2,var}}'; done | sort -k1,1 -k2,2 -T %s | datamash groupby 1,2 count 3 collapse 3 | awk -F'\\t' -v OFS='\\t' '{if ( $3 >=  %d ) {print $0}}' >> %s" % (in_tsv,working_dir,min_samples,outfile) 
	else:
		command = "for file in $(cat %s);do BASE=$(basename $file); grep -v \"^#\" $file | awk -F'\\t' -v OFS='\\t' -v var=\"$BASE\" '{if ($6 != \".\") {print $1,$2,var}}' | sed 's/^chr//g' ; done | sort -k1,1 -k2,2 -T %s | datamash groupby 1,2 count 3 collapse 3 | awk -F'\\t' -v OFS='\\t' '{if ( $3 >=  %d ) {print $0}}' >> %s" % (in_tsv,working_dir,min_samples,outfile) 

	# Submit linux command
	try:
		subprocess.run(command, shell=True)
	except subprocess.CalledProcessError as error:
		print(error)

# -------------------------
# PoN construction
# -------------------------
if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('\nTotal computing time: ' + str(round(stop - start,2)) + ' seconds')




