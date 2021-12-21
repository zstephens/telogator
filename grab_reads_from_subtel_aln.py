import os
import gzip
import argparse

def strip_polymerase_coords(rn):
	return '/'.join(rn.split('/')[:-1])

def get_sra_readname(rn):
	return rn.split(' ')[0]

def get_ccs_readname(rn):
	return rn.split(' ')[0]

def exists_and_is_nonzero(fn):
	if os.path.isfile(fn):
		if os.path.getsize(fn) > 0:
			return True
	return False

def read_fa_entry(fa_file):
	my_name = fa_file.readline().strip()
	if not my_name or my_name[0] != '>':
		return ('','')
	my_name = my_name[1:]
	my_rdat = ''
	while True:
		new_line = fa_file.readline()
		if len(new_line) and new_line[0] != '>':
			my_rdat += new_line.strip()
		else:
			fa_file.seek(fa_file.tell()-len(new_line))
			break
	return (my_name, my_rdat)

def read_fa_entry_faster(fa_file):
	my_name = fa_file.readline().strip()
	if not my_name or my_name[0] != '>':
		return ('','')
	my_name = my_name[1:]
	my_rdat = fa_file.readline().strip()
	return (my_name, my_rdat)

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='grab_subreads_from_t2t-and-subtel_aln.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('--bam',      type=str, required=True,  metavar='input.bam',        help="* Input BAM")
	parser.add_argument('--fa',       type=str, required=True,  metavar='input.fa',         help="* Raw subreads (.fa or .fa.gz)")
	parser.add_argument('--out',      type=str, required=True,  metavar='output.fa',        help="* Output reads (.fa or .fa.gz)")
	parser.add_argument('--bed',      type=str, required=True,  metavar='subtel.bed',       help="* Subtel regions")
	parser.add_argument('--samtools', type=str, required=False, metavar='samtools',         help="/path/to/samtools", default='samtools')
	parser.add_argument('--readtype', type=str, required=False, metavar='CCS / CLR / SRA',  help="Read name format",  default='CLR')
	parser.add_argument('--newlines', type=str, required=False, metavar='yes / no',         help="Read data split across multiple lines?",  default='yes')
	args = parser.parse_args()

	IN_BAM     = args.bam
	IN_READS   = args.fa
	#
	OUT_READS  = args.out
	#
	SUBTEL_BED = args.bed
	SAMTOOLS   = args.samtools
	READTYPE   = args.readtype
	#
	READ_HAS_NEWLINES = args.newlines

	if exists_and_is_nonzero(IN_BAM) == False:
		print('Error: input.bam not found.')
		exit(1)
	if exists_and_is_nonzero(IN_READS) == False:
		print('Error: input.fa not found.')
		exit(1)
	if exists_and_is_nonzero(SUBTEL_BED) == False:
		print('Error: subtel.bed not found.')
		exit(1)
	if READ_HAS_NEWLINES not in ['yes', 'no']:
		print('Error: --newlines must be yes or no')
		exit(1)

	bed_str = []
	f = open(SUBTEL_BED, 'r')
	for line in f:
		splt = line.strip().split('\t')
		bed_str.append([n for n in splt[:3]])
	f.close()

	OUT_DIR = '/'.join(OUT_READS.split('/')[:-1])
	if len(OUT_DIR) == 0:
		OUT_DIR = '.'
	OUT_DIR += '/'
	#TEMP_READNAMES = OUT_DIR + 'readnames.txt'
	TEMP_READNAMES = OUT_READS + '.readnames'

	if exists_and_is_nonzero(TEMP_READNAMES) == False:
		print('getting readnames from bam...')
		os.system('touch ' + TEMP_READNAMES)
		for i in range(len(bed_str)):
			bed_dat = bed_str[i][0] + ':' + bed_str[i][1] + '-' + bed_str[i][2]
			print('-', bed_dat)
			cmd  = SAMTOOLS + ' view ' + IN_BAM + ' ' + bed_dat + ' | cut -f 1 >> ' + TEMP_READNAMES
			os.system(cmd)

	if exists_and_is_nonzero(TEMP_READNAMES) == False:
		print('Error: failed to create readnames?')
		exit(1)

	rn_dict = {}
	f = open(TEMP_READNAMES, 'r')
	for line in f:
		if READTYPE == 'CLR':
			rn_dict[strip_polymerase_coords(line.strip())] = True
		elif READTYPE == 'SRA':
			rn_dict[get_sra_readname(line.strip())] = True
		elif READTYPE == 'CCS':
			rn_dict[get_ccs_readname(line.strip())] = True
		else:
			print('Error: unknown read type, must be: CLR or SRA')
			exit(1)
	f.close()

	if OUT_READS[-3:] == '.gz':
		f_out = gzip.open(OUT_READS, 'wt')
	else:
		f_out = open(OUT_READS, 'w')
	if IN_READS[-3:] == '.gz':
		print('getting reads from fa (gzipped)...')
		f = gzip.open(IN_READS, 'rt')
	else:
		print('getting reads from fa...')
		f = open(IN_READS, 'r')
	while True:
		if READ_HAS_NEWLINES == 'yes':
			(my_name, my_rdat) = read_fa_entry(f)
		elif READ_HAS_NEWLINES == 'no':
			(my_name, my_rdat) = read_fa_entry_faster(f)
		if not len(my_rdat):
			break
		if READTYPE == 'CLR' and strip_polymerase_coords(my_name) in rn_dict:
			f_out.write('>' + my_name + '\n' + my_rdat + '\n')
		elif READTYPE == 'SRA' and get_sra_readname(my_name) in rn_dict:
			f_out.write('>' + my_name + '\n' + my_rdat + '\n')
		elif READTYPE == 'CCS' and get_ccs_readname(my_name) in rn_dict:
			f_out.write('>' + my_name + '\n' + my_rdat + '\n')
	f.close()
	f_out.close()

	# cleanup?
	if exists_and_is_nonzero(OUT_READS) and exists_and_is_nonzero(TEMP_READNAMES):
		os.system('rm ' + TEMP_READNAMES)

if __name__ == '__main__':
	main()
