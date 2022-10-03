import argparse
import copy
import os
import random
import subprocess
import sys

# parent directory of this script
telogator_base_dir = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
sys.path.append(telogator_base_dir)

from source.tg_util import makedir, rm, exists_and_is_nonzero, read_fq_entry, index_ref, read_ref_py3
from source.tg_kmer import read_kmer_prob_matrix, sample_telomere

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23}

sorted_refs = [[n[1]+'_p', n[1]+'_q'] for n in sorted([(LEXICO_2_IND[k], k) for k in LEXICO_2_IND.keys()])]
sorted_refs = [val for sublist in sorted_refs for val in sublist]

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='simulate_reads.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-o', type=str, required=True,  metavar='output/',   help="* Output directory")
	parser.add_argument('-c', type=str, required=False, metavar='input.cfg', help="Config file", default='')
	args = parser.parse_args()

	telogator_resource_dir = telogator_base_dir + 'resources/'

	TELOGATOR_PY = telogator_base_dir + 'telogator.py'
	MERGEJOBS_PY = telogator_base_dir + 'merge_jobs.py'

	PBSIM_MODEL    = telogator_resource_dir + 'model_qc_clr'
	TELOGATOR_REF  = telogator_resource_dir + 't2t-telogator-ref.fa'
	KMER_MATRIX_FN = telogator_resource_dir + 'kmer_trans.dat'
	#
	trans_dist_dict = read_kmer_prob_matrix(KMER_MATRIX_FN, return_trans_dict=True)

	# unzip ref if necessary
	if exists_and_is_nonzero(TELOGATOR_REF+'.gz') and exists_and_is_nonzero(TELOGATOR_REF) == False:
		os.system('gunzip -c ' + TELOGATOR_REF+'.gz' + ' > ' + TELOGATOR_REF)
	# get all the reference strings we need
	ref_inds = index_ref(TELOGATOR_REF)
	ri_rev   = {ref_inds[n][0]:n for n in range(len(ref_inds))}
	ref_seqs = {}
	for ref in sorted_refs:
		my_ri = ri_rev[ref]
		ref_seqs[ref] = read_ref_py3(TELOGATOR_REF, ref_inds[my_ri])

	IN_CFG   = args.c
	OUT_BASE = args.o
	#
	makedir(OUT_BASE)

	cfg_params = {}
	#
	# default values
	#
	cfg_params['PBSIM_EXE']      = 'pbsim'
	cfg_params['PBMM2_EXE']      = 'pbmm2'
	cfg_params['SAMTOOLS']       = 'samtools'
	cfg_params['SUBTEL_CHOP']    = '10000,400000'
	cfg_params['TEL_LEN_RNG']    = '2000,10000'
	cfg_params['READLEN_SWEEP']  = '20000'
	cfg_params['ACCURACY_SWEEP'] = '0,5,10,15'
	cfg_params['COVERAGE_SWEEP'] = '20'
	cfg_params['MULTIPLICITY']   = '1'
	#
	# read in from cfg if provided
	#
	if len(IN_CFG) and exists_and_is_nonzero(IN_CFG):
		print('reading ' + IN_CFG + '...')
		f = open(IN_CFG, 'r')
		for line in f:
			line_dat = line.strip().replace(' ','')
			if len(line_dat) and '=' in line_dat:
				splt = line_dat.split('=')
				cfg_params[splt[0]] = splt[1]
		f.close()
	#
	PBSIM_EXE      = cfg_params['PBSIM_EXE']
	PBMM2_EXE      = cfg_params['PBMM2_EXE']
	SAMTOOLS       = cfg_params['SAMTOOLS']
	SUBTEL_CHOP    = [int(n) for n in cfg_params['SUBTEL_CHOP'].split(',')]
	TEL_LEN_RNG    = [int(n) for n in cfg_params['TEL_LEN_RNG'].split(',')]
	READLEN_SWEEP  = [int(n) for n in cfg_params['READLEN_SWEEP'].split(',')]
	ACCURACY_SWEEP = [int(n) for n in cfg_params['ACCURACY_SWEEP'].split(',')]
	COVERAGE_SWEEP = [int(n) for n in cfg_params['COVERAGE_SWEEP'].split(',')]
	MULTIPLICITY   = int(cfg_params['MULTIPLICITY'])
	#
	# various hardcoded simulation params
	#
	N_REF_BUFF  = 'N'*100000
	READLEN_SD  = 200
	READ_MAXDEV = 1000
	ACCURACY_SD = 0.005

	#
	# function to run PBSIM to generate simulated PacBio long reads
	#
	def run_pbsim(out_dir, ref_fa, sim_readlen, sim_accuracy, sim_coverage):
		if out_dir[-1] != '/':
			out_dir += '/'
		out_log    = out_dir + 'log.txt'
		out_prefix = 'sim'
		out_fa     = out_dir + 'sim.fa'
		PARAMS = ['--model_qc',      PBSIM_MODEL,
	              '--data-type',     'CLR',
	              '--length-mean',   str(sim_readlen),
	              '--length-sd',     str(READLEN_SD),
	              '--length-min',    str(sim_readlen - READ_MAXDEV),
	              '--length-max',    str(sim_readlen + READ_MAXDEV),
	              '--accuracy-mean', str(sim_accuracy),
	              '--accuracy-sd',   str(ACCURACY_SD),
	              '--depth',         str(sim_coverage),
	              '--prefix',        out_dir + out_prefix]
		min_read_length_n = sim_readlen - 2*READ_MAXDEV
		#
		pbsim_cmd = PBSIM_EXE + ' ' + ' '.join(PARAMS) + ' ' + ref_fa + ' > ' + out_log + ' 2>&1'
		print(pbsim_cmd)
		print('simulating reads...')
		os.system(pbsim_cmd)
		#
		f_out   = open(out_fa, 'w')
		listing = sorted([n for n in os.listdir(out_dir) if (n[-6:] == '.fastq' and n[:len(out_prefix)] == out_prefix)])
		print('consolidating reads into single fasta...')
		for i in range(len(listing)):
			f = open(out_dir + listing[i], 'r')
			while True:
				(my_r, my_q, my_n) = read_fq_entry(f)
				if my_r == '':
					break
				if 'N' in my_r:				# trim reads with Ns
					splt = my_r.split('N')
					if len(splt[0]) >= 3:
						my_r = splt[0]
					elif len(splt[-1]) >= 3:
						my_r = splt[-1]
					else:
						continue
				if len(my_r) < min_read_length_n:
					continue
				#f_out.write('>' + my_n + '\n')
				f_out.write('>' + my_n + '-' + sorted_refs[i] + '\n')
				f_out.write(my_r + '\n')
			f.close()
		f_out.close()
		#
		listing = [n for n in os.listdir(out_dir) if (n[-3:] != '.fa' and n[-4:] != '.tsv')]
		for fn in listing:
			rm(out_dir + fn)

	#
	# function to align simulated reads to the telogator reference using pbmm2
	#
	def run_pbmm2(my_reads, out_bam):
		PARAMS = ['--preset', 'SUBREAD',
		          '--sample', 'sample',
		          '--rg',     "'@RG\\tID:movie1'",
		          '-j',       '6',	# alignment threads
		          '-J',       '3',	# sorting threads
		          '--sort']
		pbmm2_cmd = PBMM2_EXE + ' align ' + TELOGATOR_REF + ' ' + my_reads + ' ' + out_bam + ' ' + ' '.join(PARAMS)
		print(pbmm2_cmd)
		print('aligning reads with pbmm2...')
		os.system(pbmm2_cmd)

	#
	# function to run telogator
	#
	def run_telogator(my_bam, out_dir, my_tlen, n_jobs=6):
		my_sam  = my_bam[:-4] + '.sam'
		sam_cmd = SAMTOOLS + ' view ' + my_bam + ' > ' + my_sam
		os.system(sam_cmd)
		sp = []
		for i in range(1,n_jobs+1):
			tel_lst = ['python3', TELOGATOR_PY, '-i', my_sam, '-o', out_dir, '--job', str(i), str(n_jobs)]
			sp.append(subprocess.Popen(tel_lst))
		exit_codes = [p.wait() for p in sp]
		print('exit codes:', exit_codes)
		com_cmd = 'python3 ' + MERGEJOBS_PY + ' -i ' + out_dir + ' -gt ' + my_tlen + ' -t p90 -rc 1 --pbsim'
		os.system(com_cmd)

	#
	# BEGIN SIMULATIONS!
	#
	for i in range(len(READLEN_SWEEP)):
		sim_ref_len = min([100000, 10*READLEN_SWEEP[i]])
		for j in range(len(ACCURACY_SWEEP)):
			for j2 in range(len(COVERAGE_SWEEP)):
				for k in range(MULTIPLICITY):
					#
					my_dir = OUT_BASE + 'sim_c-' + str(COVERAGE_SWEEP[j2]) + '_r-' + str(READLEN_SWEEP[i]) + '_a-' + str(ACCURACY_SWEEP[j]) + '_m-' + str(k) + '/'
					makedir(my_dir)
					#
					#	ARTIFICIAL SAMPLE GENOME
					#
					my_ref  = my_dir + 'ref.fa'
					my_tlen = my_dir + 'tlen.tsv'
					if exists_and_is_nonzero(my_ref) == False:
						f  = open(my_ref, 'w')
						ft = open(my_tlen, 'w')
						for ref in sorted_refs:
							ref_copy = copy.deepcopy(ref_seqs[ref])
							stc = random.randrange(SUBTEL_CHOP[0], SUBTEL_CHOP[1])
							mtl = random.randrange(TEL_LEN_RNG[0], TEL_LEN_RNG[1])
							f.write('>' + ref + '\n')
							if ref[-2:] == '_p':
								my_coord = stc
								ref_copy = ref_copy[:my_coord+sim_ref_len] + 'N'*len(ref_copy[my_coord+sim_ref_len:])
								f.write(N_REF_BUFF + sample_telomere(trans_dist_dict, mtl, tel_type='p') + ref_copy[my_coord:] + '\n')
							elif ref[-2:] == '_q':
								my_coord = len(ref_copy) - stc
								ref_copy = 'N'*len(ref_copy[:my_coord-sim_ref_len]) + ref_copy[my_coord-sim_ref_len:]
								f.write(ref_copy[:my_coord] + sample_telomere(trans_dist_dict, mtl, tel_type='q') + N_REF_BUFF + '\n')
							else:
								print('Error: check your ps and qs')
								exit(1)
							ft.write(ref + '\t' + str(my_coord) + '\t' + str(mtl) + '\n')
						ft.close()
						f.close()
					#
					#	SIMULATE READS
					#
					my_acc   = float(100 - ACCURACY_SWEEP[j])/100.
					my_reads = my_dir + 'sim.fa'
					if exists_and_is_nonzero(my_reads) == False:
						run_pbsim(my_dir, my_ref, READLEN_SWEEP[i], my_acc, COVERAGE_SWEEP[j2])
					#
					#	PBMM2
					#
					my_bam = my_dir + 'aln.bam'
					if exists_and_is_nonzero(my_bam) == False:
						run_pbmm2(my_reads, my_bam)
					#
					#	TELOGATOR
					#
					tg_dir = my_dir + 'telogator/'
					makedir(tg_dir)
					my_results = tg_dir + 'results.tsv'
					if exists_and_is_nonzero(my_results) == False:
						run_telogator(my_bam, tg_dir, my_tlen)

if __name__ == '__main__':
	main()
