import timeit
import collections
import argparse
import glob
import os
import time

def updatefields(REF,INFO,BC,chrs,pos,i,line): #It should check that if the line is null it should just update some flag
	line = line.split('\t')
	CHROM, POS, ref, info, bc = line

	if line != '':

		if CHROM==chrs[i] and  pos[i]>=int(POS):
			pass
		else:
			chrs[i] = CHROM
			pos[i] = int(POS)
			REF[i] = ref
			INFO[i] = info
			BC[i] = bc
		
	return (REF,INFO,BC,chrs,pos)

def should_we_go(lines):
	#False when all lines are emtpy
	for i in lines:
		if i != "":
			return True
	return False    


def read_one_line(chr_file,pos_file,cur_chr,cur_pos):
	#It checks that the current file is correct to read one more line
	#meanstaht the chr should be the same and the position should be <= to the current position
	if chr_file == cur_chr and pos_file <= cur_pos:
		return True
	else :
		return False

def check_chr(cur_chr,chrs):
	#it returns true if at least one chr is equal to the current chr
	for c in chrs:
		if c == cur_chr:
			return True
	return False

def sort_set(LIST):
	while ('NA' in LIST):
		LIST.remove('NA')
	
	while (len(LIST) > 1 and '.' in LIST):
		LIST.remove('.')
	
	counter=collections.Counter(LIST)
	sort_set_list = [pair[0] for pair in sorted(counter.items(), key=lambda item: item[1], reverse=True)]
	return('|'.join(sort_set_list))

def write_me(chrs,pos, REF, INFO, BC,cur_chr,outfile,low_pos):
	# extract the correct field to write to the output file
	# for each of the files check those at the lowest position possible
	# and produce the output from them

	low_pos_REF = list()

	low_pos_INFO = list()

	# Where we save all the information per cell type
	low_pos_CellType = list()

	for i in range(len(chrs)):
		if chrs[i]==cur_chr and pos[i]==low_pos:
			low_pos_REF.append(str(REF[i]))
			low_pos_INFO.append(str(INFO[i]))
			low_pos_CellType.append(str(BC[i]))
		else:
			low_pos_REF.append('NA')
			low_pos_INFO.append(str('NA'))
			low_pos_CellType.append('NA')
		
	info = '\t'.join([str(cur_chr),str(low_pos),str(low_pos), sort_set(low_pos_REF), sort_set(INFO)])
	
	out_cell_types = info + '\t' + '\t'.join(low_pos_CellType)
	outfile.write(out_cell_types + '\n')
	
	return None


def calc_low(chrs,pos,cur_chr,cur_pos): #calculate the lowest current position
	
	pos_in_cur_chr= [pos[x] for x in range(len(chrs)) if chrs[x] == cur_chr  and pos[x]>cur_pos]
	
	try:
		low_pos = min(pos_in_cur_chr)
		if low_pos <= cur_pos:
			print ('error')
			print (low_pos)
			raise
	except:
		print ('error')
		print (chrs)
		print (pos)
		print (cur_chr)
		print (cur_pos)
		raise
	return  low_pos

def update_chr(chrs):
	# returns the current chr from the  reads of the files
	# update the chromosome to move across files
	new_chrs = sorted([s for s in chrs if s])
	new_chr = new_chrs[0]
	return new_chr


def merge_cell_types_files(infiles,outfile):

	fileids = list()
	chrs = list()
	pos = list()
	lines=list()
	
	# Reference nucleotide
	REF = list()
	# Coverage
	INFO = list()
	# Base counts
	BC = list()


	header_lines = 9 # Number of header lines of the input tsvs

	# Common header for all files
	header=['#CHROM','Start','End','REF','INFO']
	date = time.strftime("%d/%m/%Y")## dd/mm/yyyy format
	DATE="##fileDate=%s" % date
	CONCEPTS="""##INFO=DP,Description="Depth of coverage">\n##INFO=NC,Description="Number of different cells">\n##INFO=CC,Description="Cell counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BC,Description="Base counts [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BQ,Description="Base quality sums [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BCf,Description="Base counts in forward reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">\n##INFO=BCr,Description="Base counts in reverse reads [A:C:T:G:I:D:N:O], where D means deletion, I insertion and O other type of character">"""

	for id_ in range(len(infiles)):
		
			tmp=open(infiles[id_],'r')
			fileids.append(tmp)
			for i in range(header_lines):
				useless = fileids[id_].readline() #skip the header
			useless = fileids[id_].readline()
			lines.append(useless)
			curline=useless.strip()
			
			# Create empty sites
			REF.append(0)
			INFO.append(0)
			BC.append(0)

			chrs.append('x')
			pos.append(0)
			
			# Get cell types and append them to header
			CELL_TYPE = os.path.basename(infiles[id_]).split('.')[-2]
			header.append(CELL_TYPE)

			updatefields(REF,INFO,BC,chrs,pos,id_,curline)

	cur_chr=1
	cur_pos=0

	# Create output file
	outfile= open(outfile,'w')
	
	# Write header
	outfile.write(DATE + '\n')
	outfile.write(CONCEPTS + '\n')
	outfile.write('\t'.join(header)+'\n')

	# Check if we should countinue processing files
	go = should_we_go(lines)
	
	while go:
		for i in range(len(infiles)):
		
			while read_one_line(chrs[i],pos[i],cur_chr,cur_pos):
				lines[i]=fileids[i].readline().strip()
				
				if lines[i]=="":
					chrs[i] = ''
					pos[i] = -1
					break
		
				REF,INFO,BC,chrs,pos = updatefields(REF,INFO,BC,chrs,pos,i,lines[i]) # updates the fields for the current line

		go = should_we_go(lines)
		if check_chr(cur_chr,chrs):
			
			if go:
				cur_pos = calc_low(chrs,pos,cur_chr,cur_pos) #calculate the lowest current position
				write_me(chrs,pos, REF, INFO, BC, cur_chr, outfile, cur_pos)
		else:
			#print lines
			if should_we_go(lines):
			  cur_chr = update_chr(chrs) #routine to update the current chr, should check that's the same for all
			  cur_pos = 0
			
		pass        

	outfile.close()

def initialize_parser():
	parser = argparse.ArgumentParser(description='Script to merge the cell/base counts tsv files per cell type in only one')
	parser.add_argument('--tsv_folder', type=str, default=1, help='Path to the directory containing the base count files in tsv format for each cell type. All tsv files in the directory will be used. Avoid not desired tsv files in this folder', required = True)   
	parser.add_argument('--outfile', help='Output file name', required = True)
	return (parser)

def main():

	# 1. Arguments
	parser = initialize_parser()
	args = parser.parse_args()

	tsv_folder = args.tsv_folder
	out_file = args.outfile

	# 1. Merge cell types files and select potential het sites for phasing
	print ('-----------------------------------------------------------')
	print ('1. Merging cell types in a unique tsv file')
	print ('-----------------------------------------------------------\n')

	# Checking the number of tsv files in the folder
	infiles = glob.glob(tsv_folder + '/*.tsv')

	if len(infiles) < 1:
		raise RuntimeError('No tsv files found')
	else:
		print (str(len(infiles)) + ' tsv files found\n')

	# Merging files
	merge_cell_types_files(infiles,out_file)

	print ('Done...\n')


if __name__ == '__main__':
	start = timeit.default_timer()
	main()
	stop = timeit.default_timer()
	print ('Time: ' + str(round(stop - start,2)) + ' seconds')

