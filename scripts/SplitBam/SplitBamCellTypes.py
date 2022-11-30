import pysam
import pandas as pd
import argparse
import timeit
import sys
import numpy as np

def meta_to_dict(txt,tissue):
	metadata = pd.read_csv(txt, delimiter = "\t")

	# Clean index column
	metadata['Index_clean'] = metadata['Index'].str.replace('-.*$','',regex=True)

	# If tissue provided, append tissue id to cell type to be recognised in downstream analysis
	if tissue == None:
		metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ','_',regex=True)
	else:
		tissue = tissue.replace(" ", "_")
		metadata['Cell_type_clean'] = metadata['Cell_type'].str.replace(' ','_',regex=True)
		metadata['Cell_type_clean'] = str(tissue) + '__' + metadata['Cell_type_clean'].astype(str)

	# Create dicitionary with cell types and cell barcodes
	DICT = metadata.set_index('Index_clean')['Cell_type_clean'].to_dict()
	ALL_CELL_TYPES = metadata['Cell_type_clean'].unique()

	del metadata

	return(DICT, ALL_CELL_TYPES)


def split_bam(bam, txt, outdir,donor,tissue,max_NM,max_NH,min_MAPQ,n_trim):

	start = timeit.default_timer()

	# 1. Transform table to dictionary
	DICT, ALL_CELL_TYPES = meta_to_dict(txt,tissue)

	if len(DICT) < 1:
		print ('Warning: No cell barcodes found in the --meta file')
		sys.exit()

	# 2. Open infile
	infile=pysam.Samfile(bam, "rb") 

	# 3. Create and open out bam files
	DICT_files = {}
	for cell_type in ALL_CELL_TYPES:
		outfile="{}/{}.{}.bam".format(outdir, donor, cell_type)
		outfile_wb = pysam.AlignmentFile(outfile, "wb",template=infile)
		DICT_files[cell_type] = outfile_wb

	# 4. Start read counts
	total_reads = 0
	FILTER_dict = {'Total_reads' : 0,'Pass_reads' : 0, 'CB_not_found' : 0, 'CB_not_matched' : 0} # To store filter reasons
	
	# 5. Check reads and split them in bam files
	for read in infile.fetch():
		FILTER_dict['Total_reads'] += 1 # To count the total number of reads analysed
		total_reads += 1

		if (total_reads % 5000000 == 0):
			print ('Number of reads already processed: ' + str(total_reads))

		# Check if CB tag is present
		try:
			barcode = read.opt("CB")
		except:
			FILTER2 = 'CB_not_found'
			FILTER_dict[FILTER2] += 1
			continue 


		# Check if CB code matches with an annotated cell type
		barcode = barcode.split("-")[0]
		try:
			CELL_TYPE = DICT[barcode]
		except:
			FILTER2 = 'CB_not_matched'
			FILTER_dict[FILTER2] += 1
			continue 

		# Final filters
		FILTER = [] 
		# Check number of mismatches in the read
		if (max_NM != None):
			try:
				if (read.opt("nM") > max_NM):
					FILTER.append('nM')
			except KeyError as e:
				FILTER.append('nM_not_found')
		# Check number of hits
		if (max_NH != None):
			try:
				if (read.opt("NH") > max_NH):
					FILTER.append('NH')
			except KeyError as e:
				FILTER.append('NH_not_found')

		# Check mapping quality
		if (min_MAPQ > 0):
			try:
				if (read.mapq < min_MAPQ):
					FILTER.append('MAPQ')
			except:
				FILTER.append('MAPQ_not_found')

		# Making a decision about the read (filtered or not)
		if (len(FILTER) > 0): # If there are reasons to filter read
			FILTER2 = ';'.join(FILTER)
			
			FILTER_dict[FILTER2] = FILTER_dict.get(FILTER2, 0) + 1
			continue
		else:
			FILTER_dict['Pass_reads'] += 1

		# Only for PASS reads
		# Trim last and first bases of the read (reduce quality) if specified
		# It does not consider the soft-clip bases at the beginning/end of the read for the counting
		# It also considers the expected adapter lengths (up to 30) of 10x library prep to remove bases in long soft-clip sequences
		if (n_trim > 0):
			CIGAR_tuple = read.cigartuples
			if (CIGAR_tuple != None and len(CIGAR_tuple) > 1):
				# Check the number of bases to trim at the beginning of the read
				# Extension of the soft-clipped bases
				if (CIGAR_tuple[0][0] == 4):
					# As there are some library preparation adapters inside the 10x prep protocol that 
					# are not properly removed (up to 30 bps), if we observed that the softclip bases are around more than 20, we 
					# prefer to be conservative and remove the expected 30 pbs + the n_trim parameter
					if (CIGAR_tuple[0][1] >= 20 and CIGAR_tuple[0][1] < 30):
						trim_start = 30 + n_trim
					else:
						trim_start = CIGAR_tuple[0][1] + n_trim
				else:
					trim_start = n_trim

				# Check the number of bases to trim at the end of the read
				# Extension of the soft-clipped bases
				if (CIGAR_tuple[-1][0] == 4):	
					# As there are some library preparation adapters inside the 10x prep protocol that 
					# are not properly removed (up to 30 bps), if we observed that the softclip bases are around more than 20, we 
					# prefer to be conservative and remove the expected 30 pbs + the n_trim parameter
					if (CIGAR_tuple[-1][1] >= 20 and CIGAR_tuple[-1][1] < 30):
						trim_end = 30 + n_trim
					else:
						trim_end = CIGAR_tuple[-1][1] + n_trim
				else:
					trim_end = n_trim
			else:
				trim_start = n_trim
				trim_end = n_trim

			# Set base qualities to 0 at the beginning and the end of the read if specified
			# Get first and last qualities
			Q_indexes = list(np.arange(trim_start)) + list((np.arange(trim_end) + 1)*-1)
			Q = read.query_qualities
			for Qi in Q_indexes:
				Q[Qi] = 0

			# Substitute qualities
			read.query_qualities = Q

		# Print passed read
		DICT_files[CELL_TYPE].write(read)

	# 6. Close opened files
	infile.close()
	for key in DICT_files.keys():
		DICT_files[key].close()

	# 7. Get and create report
	outfile_report = outfile="{}/{}.report.txt".format(outdir, donor, cell_type)
	stop = timeit.default_timer()
	endtime = round((stop - start),2)
	FILTER_dict['Total_time'] = endtime

	data_df = pd.DataFrame([FILTER_dict])
	data_df.to_csv(outfile_report, index = False, sep = '\t')

	# 8. Index bam files
	for cell_type in ALL_CELL_TYPES:
		bam="{}/{}.{}.bam".format(outdir, donor, cell_type)
		pysam.index(bam)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Split alignment file into cell type specific BAMs')
	parser.add_argument('--bam', type=str, default=1, help='BAM file to be analysed (Sorted by coordinate)', required = True)
	parser.add_argument('--meta', type=str, default=1, help='Metadata file mapping cell barcodes to cell type information', required = True)
	parser.add_argument('--id', type=str, default = 'Sample', help='Sample ID', required = False)
	parser.add_argument('--max_nM', type=int, default = None, help='Maximum number of mismatches permitted to consider reads for analysis. By default, this filter is switched off, although we recommed using --max_nM 5. If applied, this filter requires having the nM tag in the bam file. [Default: Switched off]', required = False)
	parser.add_argument('--max_NH', type=int, default = None, help='Maximum number of alignment hits permitted to consider reads for analysis. By default, this filter is switched off, although we recommend using --max_NH 1. This filter requires having the NH tag in the bam file. [Default: Switched off]', required = False)
	parser.add_argument('--min_MQ', type=int, default = 255, help='Minimum mapping quality required to consider reads for analysis. Set this value to 0 to switch this filter off. --min_MQ 255 is recommended for RNA data, and --min_MQ 30 for DNA data. [Default: 255]', required = False)
	parser.add_argument('--n_trim', type=int, default = 0, help='Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read [Default: 0]', required = False)
	parser.add_argument('--outdir', default = '.', help='Out directory', required = False)
	return (parser)


def main():
	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	bam = args.bam
	outdir = args.outdir
	txt = args.meta
	donor = args.id
	tissue = None
	max_NM = args.max_nM
	min_MAPQ = args.min_MQ
	n_trim = args.n_trim
	max_NH = args.max_NH

	# 2. Split bam file
	split_bam(bam, txt, outdir,donor,tissue,max_NM,max_NH,min_MAPQ,n_trim)


#-------------
# Execute code
#-------------

if __name__ == '__main__':
	main()


