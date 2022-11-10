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


def split_bam(bam, txt, outdir,donor,tissue,max_NM,min_MAPQ,n_trim):

	# max_NM = 1000
	# min_MAPQ = 255
	# n_trim = 3

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
	not_CB = 0
	not_found = 0
	pass_reads = 0
	parameter_filter = 0
	parameter_filter_not_found = 0

	# 5. Check reads and split them in bam files
	for read in infile.fetch():
		total_reads = total_reads + 1

		if (total_reads % 5000000 == 0):
			print ('Number of reads already processed: ' + str(total_reads))

		try:
			barcode = read.opt("CB")
		except:
			not_CB = not_CB + 1
			continue

		barcode = barcode.split("-")[0]

		try:
			CELL_TYPE = DICT[barcode]
		except:
			not_found = not_found + 1
			continue

		# Final filters
		if (max_NM < 1000 or min_MAPQ > 0):
			try:
				if (read.opt("nM") <= max_NM and read.mapq >= min_MAPQ and read.opt("NH") < 2):
					pass_reads = pass_reads + 1
				else:
					parameter_filter = parameter_filter + 1
					continue
			except:
				parameter_filter_not_found = parameter_filter_not_found + 1
				continue
		else:
			pass_reads = pass_reads + 1

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

	# 7. Get report
	outfile_report = outfile="{}/{}.report.txt".format(outdir, donor, cell_type)
	stop = timeit.default_timer()
	endtime = round((stop - start),2)
	data = [{'Total_reads': total_reads, 'Pass_reads': pass_reads, 'Reads_without_cell_type': not_found, 'Reads_without_CB': not_CB, 'Total_time':endtime}]
	data_df = pd.DataFrame.from_records(data)

	data_df.to_csv(outfile_report, index = False, sep = '\t')

	# 8. Index bam files
	for cell_type in ALL_CELL_TYPES:
		bam="{}/{}.{}.bam".format(outdir, donor, cell_type)
		pysam.index(bam)

def initialize_parser():
	parser = argparse.ArgumentParser(description='Split alignment file into cell type specific BAMs')
	parser.add_argument('--bam', type=str, default=1, help='BAM file to be analysed (Sorted by coordinates)', required = True)
	parser.add_argument('--meta', type=str, default=1, help='Metadata file mapping cell barcodes to cell type information', required = True)
	parser.add_argument('--id', type=str, default = 'Sample', help='Sample ID', required = False)
	parser.add_argument('--tissue', type=str, default = None, help='Tissue name. Recommended when different tissues from the same individual are analysed', required = False)
	parser.add_argument('--max_nM', type=int, default = 1000, help='Maximum number of mismatches permitted to consider a read for further analysis [Default: 1000]', required = False)
	parser.add_argument('--min_MQ', type=int, default = 255, help='Minimum mapping quality required to consider a read for further analysis [Default: 255]', required = False)
	parser.add_argument('--n_trim', type=int, default = 0, help='Number of bases trimmed by setting the base quality to 0 at the beginning and end of the read [Default: 0]', required = False)
	parser.add_argument('--outdir', default = '.', help='Output directory', required = False)
	return (parser)


def main():
	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	bam = args.bam
	outdir = args.outdir
	txt = args.meta
	donor = args.id
	tissue = args.tissue
	max_NM = args.max_nM
	min_MAPQ = args.min_MQ
	n_trim = args.n_trim

	# 2. Split bam file
	split_bam(bam, txt, outdir,donor,tissue,max_NM,min_MAPQ,n_trim)


#-------------
# Execute code
#-------------

if __name__ == '__main__':
	main()


