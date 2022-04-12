import argparse
import copy
import gzip
import pathlib
import multiprocessing

import numpy as np

from source.tg_kmer   import get_telomere_kmer_density, get_telomere_regions
from source.tg_reader import TG_Reader
from source.tg_util   import posmax

TEL_WINDOW_SIZE = 100
MIN_TEL_SCORE   = 100
BTWN_TEL_SCORE  = 0.8
MAX_EDGE_NONTEL = 1000
#

def process_reads(MYJOB, NJOBS, IN_READS, KMER_DICT, P_VS_Q_AMP_THRESH, MODE, return_dict):
	#
	nreads_processed   = 0
	return_dict[MYJOB] = []
	#
	for job_i in range(len(IN_READS)):
		rdat_list = IN_READS[job_i]

		all_tel_regions = []
		for (rdat, offset) in rdat_list:
			# compute telomere kmer density
			(td_p_e0, td_p_e1) = get_telomere_kmer_density(rdat, KMER_DICT['p'], TEL_WINDOW_SIZE)
			(td_q_e0, td_q_e1) = get_telomere_kmer_density(rdat, KMER_DICT['q'], TEL_WINDOW_SIZE)
			# estimate telomere boundaries
			(p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
			#
			all_tel_regions.extend([n[0]+offset, n[1]+offset, n[2]] for n in tel_regions)
		tel_regions = all_tel_regions

		#
		# just want any tel sequence, don't care if read terminates in it or not
		#
		if MODE == 'unanchored':
			my_tel_regions = [tuple(n) for n in all_tel_regions if n[2] != None]
			for n in my_tel_regions:
				return_dict[MYJOB] += [(job_i, n[0], n[1], n[2])]

		#
		# normal telogator logic (only consider telomeres found at the end of reads)
		#
		elif MODE == 'flank':
			score_scalars = np.ones(len(tel_regions))
			for i in range(1, len(tel_regions)-1):
				if tel_regions[i][2] == None and tel_regions[i-1][2] == tel_regions[i+1][2]:
					score_scalars[i] = BTWN_TEL_SCORE

			#
			my_tel_len_p  = 0
			my_tel_len_q  = 0
			my_plot_title = ''
			my_tel_aln    = []
			my_tel_len    = None
			my_tel_conc   = [{'p':0, 'q':0, None:0}, {'p':0, 'q':0, None:0}]
			#
			#	get tel lengths, from the left (p-arm)
			#
			if tel_regions[0][2] == None and tel_regions[0][1] - tel_regions[0][0] > MAX_EDGE_NONTEL:
				pass
			else:
				my_score = []
				for i in range(len(tel_regions)):
					my_s = (tel_regions[i][1] - tel_regions[i][0]) * score_scalars[i]
					if tel_regions[i][2] == 'p':
						my_score.append(my_s)
					elif tel_regions[i][2] == 'q':	# not too sure how to handle this
						my_score.append(my_s)		# should this be -my_s ??
					elif tel_regions[i][2] == None:
						my_score.append(-my_s)
				cum_score = np.cumsum(my_score)
				max_i     = posmax(cum_score)
				if cum_score[max_i] >= MIN_TEL_SCORE:
					my_tel_len_p   = int(tel_regions[max_i][1] + TEL_WINDOW_SIZE/2)
					#print('P-ARM TEL:', my_tel_len_p, [int(n) for n in cum_score.tolist()], max_i, '\n')
					for i in range(0, max_i+1):
						#print('tel_regions['+str(i)+'] =', tel_regions[i])
						my_tel_conc[0][tel_regions[i][2]] += abs(tel_regions[i][1] - tel_regions[i][0])
			#
			#	get tel lengths, from the right (q-arm)
			#
			if tel_regions[-1][2] == None and tel_regions[-1][1] - tel_regions[-1][0] > MAX_EDGE_NONTEL:
				pass
			else:
				my_score = []
				for i in range(len(tel_regions)):
					my_s = (tel_regions[i][1] - tel_regions[i][0]) * score_scalars[i]
					if tel_regions[i][2] == 'q':
						my_score.append(my_s)
					elif tel_regions[i][2] == 'p':	# not too sure how to handle this
						my_score.append(my_s)		# should this be -my_s ??
					elif tel_regions[i][2] == None:
						my_score.append(-my_s)
				cum_score = np.cumsum(my_score[::-1])[::-1]
				max_i     = posmax(cum_score)
				if cum_score[max_i] >= MIN_TEL_SCORE:
					my_tel_len_q   = int(tel_regions[-1][1] - tel_regions[max_i][0] + TEL_WINDOW_SIZE/2)
					#print('Q-ARM TEL:', my_tel_len_q, [int(n) for n in cum_score.tolist()], max_i, '\n')
					for i in range(max_i, len(tel_regions)):
						#print('tel_regions['+str(i)+'] =', tel_regions[i])
						my_tel_conc[1][tel_regions[i][2]] += abs(tel_regions[i][1] - tel_regions[i][0])

			#
			# sort out which tel to pick
			#
			which_to_choose = None
			if my_tel_len_p > 0 and my_tel_len_q > 0:
				tel_frac_p =  my_tel_len_p / float(len(rdat))
				tel_frac_q =  my_tel_len_q / float(len(rdat))
				# try to recover anchors in cases where we cruised through short subtel sequence
				if tel_frac_p > MAXIMUM_TEL_FRAC and tel_frac_q <= MAXIMUM_TEL_FRAC:
					which_to_choose = 'q'
				elif tel_frac_q > MAXIMUM_TEL_FRAC and tel_frac_p <= MAXIMUM_TEL_FRAC:
					which_to_choose = 'p'
				# otherwise pick whichever one is longest
				elif my_tel_len_p >= my_tel_len_q:
					which_to_choose = 'p'
				else:
					which_to_choose = 'q'
			else:
				if my_tel_len_p >= my_tel_len_q:
					which_to_choose = 'p'
				else:
					which_to_choose = 'q'
			#
			if which_to_choose == 'p':
				my_tel_len_q  = 0
				my_tel_len    = my_tel_len_p
				my_plot_title = ' : p-arm tel = ' + str(my_tel_len_p)
				my_tel_aln    = [[0, my_tel_len_p, 'tel-p', None, None, 'FWD']]
				my_tel_conc   = my_tel_conc[0]
			elif which_to_choose == 'q':
				my_tel_len_p  = 0
				my_tel_len    = my_tel_len_q
				my_plot_title = ' : q-arm tel = ' + str(my_tel_len_q)
				my_tel_aln    = [[len(rdat)-my_tel_len_q, len(rdat), 'tel-q', None, None, 'REV']]
				my_tel_conc   = my_tel_conc[1]

			if my_tel_len > 0:
				out_readname = 'tlen-' + which_to_choose + '-' + str(my_tel_len)
				return_dict[MYJOB] += [(job_i, out_readname)]
				print('---', out_readname)

		nreads_processed += 1
		print(MYJOB, job_i+1, '/', len(IN_READS))
		#if nreads_processed >= 30:
		#	break

#
#
#
def main(raw_args=None):
	parser = argparse.ArgumentParser(description='search_reads_for_tel.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',          type=str,   required=True,  metavar='in.fa',             help="* Input reads")
	parser.add_argument('-o',          type=str,   required=True,  metavar='out.fa',            help="* Output reads containing telomeres")
	parser.add_argument('-k',          type=str,   required=False, metavar='default_kmers.tsv', help="Telomere kmers to search for",           default='')
	parser.add_argument('-l',          type=int,   required=False, metavar='5000',              help="Min read length",                        default=5000)
	parser.add_argument('-s',          type=int,   required=False, metavar='0',                 help="Search this many bp into either side of the read (0 = entire read)", default=0)
	parser.add_argument('-p',          type=float, required=False, metavar='0.5',               help="Telomere signal threshold (0-1)",        default=0.5)
	parser.add_argument('--processes', type=int,   required=False, metavar='6',                 help="Number of multiprocesses to use",        default=6)
	parser.add_argument('--mode',      type=str,   required=False, metavar='flank',             help="flank / unanchored",                     default='flank')
	args = parser.parse_args()

	INPUT_FASTA       = args.i
	OUTPUT_FASTA      = args.o
	P_VS_Q_AMP_THRESH = args.p
	MINIMUM_READ_LEN  = args.l
	READLEN_BUFF      = args.s
	NUM_PROCESSES     = args.processes

	MODE = args.mode
	if MODE not in ['flank', 'unanchored']:
		print('Error: --mode must be flank or unanchored')
		exit(1)

	#
	KMER_FILE = args.k
	KMER_DICT = {'p':[], 'q':[]}
	if KMER_FILE == '':
		print('using default telomere kmers.')
		sim_path = pathlib.Path(__file__).resolve().parent
		kmer_fn  = str(sim_path) + '/resources/default_kmers.tsv'
		f = open(kmer_fn,'r')
		for line in f:
			splt = line.strip().split('\t')
			KMER_DICT[splt[0].lower()].append(splt[1])
		f.close()
	else:
		fn_suffix = KMER_FILE.split('/')[-1]
		print('using user-specified kmer list:', fn_suffix)
		if exists_and_is_nonzero(KMER_FILE):
			f = open(KMER_FILE,'r')
			for line in f:
				splt = line.strip().split('\t')
				KMER_DICT[splt[0].lower()].append(splt[1])
			f.close()
		else:
			print('Error: kmer list not found')
			exit(1)

	#
	my_reader = TG_Reader(INPUT_FASTA)
	#
	all_readdat = []
	nreads_read = 0
	#
	while True:
		read_dat = my_reader.get_next_read()
		if not read_dat[0]:
			break
		if len(read_dat[1]) >= MINIMUM_READ_LEN:
			all_readdat.append((read_dat[0], read_dat[1]))
			nreads_read += 1
	print('found', nreads_read, 'reads (>=' + str(MINIMUM_READ_LEN) + 'bp).')
	my_reader.close()

	#
	#
	#
	chunk_size  = len(all_readdat)//NUM_PROCESSES + 1
	manager     = multiprocessing.Manager()
	return_dict = manager.dict()
	processes   = []
	for i in range(NUM_PROCESSES):
		job_rdat = []
		for j in range(i*chunk_size, min([(i+1)*chunk_size, len(all_readdat)])):
			my_rdat = all_readdat[j][1]
			if READLEN_BUFF > 0 and len(my_rdat) > 2*READLEN_BUFF:
				job_rdat.append([(my_rdat[:READLEN_BUFF], 0), (my_rdat[-READLEN_BUFF:], len(my_rdat)-READLEN_BUFF)])
			else:
				job_rdat.append([(my_rdat, 0)])
		p = multiprocessing.Process(target=process_reads, args=(i+1, NUM_PROCESSES, job_rdat, KMER_DICT, P_VS_Q_AMP_THRESH, MODE, return_dict))
		processes.append(p)
	for i in range(NUM_PROCESSES):
		processes[i].start()
	for i in range(NUM_PROCESSES):
		processes[i].join()

	#print('RETURN:', return_dict)

	#
	#
	#
	out_dat = []
	for i in range(NUM_PROCESSES):
		offset = i*chunk_size
		#
		if MODE == 'flank':
			for (read_i, tel_readname) in return_dict[i+1]:
				my_i = offset + read_i
				print(i+1, my_i, all_readdat[my_i][0], '-->', tel_readname, len(all_readdat[my_i][1]))
				out_dat.append('>' + tel_readname + '_' + all_readdat[my_i][0] + '\n')
				out_dat.append(all_readdat[my_i][1] + '\n')
		#
		elif MODE == 'unanchored':
			for (read_i, tel_start, tel_end, tel_type) in return_dict[i+1]:
				my_i = offset + read_i
				out_dat.append('>' + all_readdat[my_i][0] + '_' + str(tel_start) + '_' + str(tel_end) + '_' + str(tel_type) + '\n')
				out_dat.append(all_readdat[my_i][1][tel_start:tel_end] + '\n')
	#
	#
	#
	print('writing output fasta...')
	if OUTPUT_FASTA[-3:] == '.gz':
		f_out = gzip.open(OUTPUT_FASTA, 'wt')
	else:
		f_out = open(OUTPUT_FASTA, 'w')
	for i in range(len(out_dat)):
		for n in out_dat[i]:
			f_out.write(n)
	f_out.close()

if __name__ == '__main__':
	main()
