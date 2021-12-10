import argparse
import copy
import pickle
import sys

import numpy as np

from source.tg_kmer import get_telomere_kmer_density, get_telomere_regions
from source.tg_plot import get_read_alignment_polygons, plot_all_read_data
from source.tg_util import exists_and_is_nonzero, makedir, parse_cigar, parse_read, repeated_matches_trimming, RC, posmax

#
TEL_WINDOW_SIZE  = 100
MIN_TEL_SCORE    = 100
BTWN_TEL_SCORE   = 0.8
MAX_EDGE_NONTEL  = 1000
#
# toss reads if too much of the non-telomere sequence couldn't be aligned anywhere
MAXIMUM_UNEXPLAINED_FRAC = 0.7
#
# toss reads if the median |p_vs_q| of the non-telomere sequence is above this value
MAX_NONTEL_MEDIAN_KMER_DENSITY = 0.25

#
# ANCHORING_STRATEGY = 'largest'     - anchor tels onto largest non-tel alignment
# ANCHORING_STRATEGY = 'closest'     - anchor tels onto nearest non-tel alignment
#
# MATCH_TRIM_STRATEGY = 'largest'    - prioritize largest alignments when trimming overlaps
# MATCH_TRIM_STRATEGY = 'mapq'       - prioritize alignments with highest MAPQ when trimming overlaps
# MATCH_TRIM_STRATEGY = 'none'       - do not trim overlaps (use at your own risk for mm2 BAMs)
#

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='Telogator v1.0', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',    type=str,   required=True,  metavar='input.sam (or merged_aln.p or - for stdin)', help="* Long reads aligned to subtel ref (- for stdin)")
	parser.add_argument('-o',    type=str,   required=True,  metavar='output/',                                    help="* Path to output directory")
	parser.add_argument('-l',    type=int,   required=False, metavar='[5000]',            default=5000,            help="Minimum read length")
	parser.add_argument('-t',    type=float, required=False, metavar='[0.9]',             default=0.9,             help="Maximum fraction of read that can be tel")
	parser.add_argument('-p',    type=float, required=False, metavar='[0.5]',             default=0.5,             help="Telomere signal threshold (0-1)")
	#
	parser.add_argument('--sa',  type=str,   required=False, metavar='[largest]',         default='largest',       help="Subtel/tel anchoring strategy")
	parser.add_argument('--sm',  type=str,   required=False, metavar='[mapq]',            default='mapq',          help="Repeated matches trimming strategy")
	#
	parser.add_argument('--job', type=int,   required=False, metavar=('my_job','n_jobs'), default=(0,0),           help='Job splitting for parallelization', nargs=2)
	parser.add_argument('--plot',            required=False, action='store_true',         default=False,           help='Create read plots')
	parser.add_argument('--debug',           required=False, action='store_true',         default=False,           help='Print results for each read as its processed')
	args = parser.parse_args()

	INPUT_SAM  = args.i
	INPUT_TYPE = None
	if INPUT_SAM == '-':
		INPUT_TYPE = 'stdin'
	else:
		if exists_and_is_nonzero(INPUT_SAM) == False:
			print('Error reading -i input file')
			exit(1)
		if INPUT_SAM[-4:].lower() == '.sam':
			INPUT_TYPE = 'sam'
		elif INPUT_SAM[-2:].lower() == '.p':
			INPUT_TYPE = 'pickle'
		else:
			print('Error: unknown -i file type')
			exit(1)

	# example for splitting into 2 jobs:
	#    --job 1 2
	#    --job 2 2
	(MYJOB, NJOBS) = args.job
	if MYJOB == 0:
		MYJOB = 1
		NJOBS = 1
	else:
		print('running job', MYJOB, 'of', NJOBS)

	OUT_PLOT_DIR = args.o
	if OUT_PLOT_DIR[-1] != '/':
		OUT_PLOT_DIR += '/'
	makedir(OUT_PLOT_DIR)
	OUT_PICKLE = OUT_PLOT_DIR + 'tel-data' + '_job-'+str(MYJOB)+'-'+str(NJOBS) + '.p'

	#
	P_VS_Q_AMP_THRESH = args.p
	MINIMUM_READ_LEN  = args.l
	MAXIMUM_TEL_FRAC  = args.t
	#
	ANCHORING_STRATEGY  = args.sa
	MATCH_TRIM_STRATEGY = args.sm

	PLOT_READS  = args.plot
	PRINT_DEBUG = args.debug

	#
	#
	#
	# [readpos_start, readpos_end, ref, pos_start, pos_end, orientation, mapq, read_sequence]
	#
	ALIGNMENTS_BY_RNAME = {}
	sorted_keys         = {}
	#
	if INPUT_TYPE == 'stdin':
		print('reading stdin input...')
		if not sys.stdin.isatty():
			input_stream = sys.stdin
		else:
			print('Error: did not receive stdin input')
			exit(1)
		for line in input_stream:
			if len(line):
				#
				[rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(line)
				#
				if rnm not in sorted_keys:
					sorted_keys[rnm] = (ref_key, pos)
				elif ref_key < sorted_keys[rnm][0] or (ref_key == sorted_keys[rnm][0] and pos < sorted_keys[rnm][1]):
					sorted_keys[rnm] = (ref_key, pos)
				if rnm not in ALIGNMENTS_BY_RNAME:
					ALIGNMENTS_BY_RNAME[rnm] = []
				ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
	#
	elif INPUT_TYPE == 'sam':
		print('reading SAM input...')
		f = open(INPUT_SAM, 'r')
		for line in f:
			#
			[rnm, ref_key, pos, read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat] = parse_read(line)
			#
			if rnm not in sorted_keys:
				sorted_keys[rnm] = (ref_key, pos)
			elif ref_key < sorted_keys[rnm][0] or (ref_key == sorted_keys[rnm][0] and pos < sorted_keys[rnm][1]):
				sorted_keys[rnm] = (ref_key, pos)
			if rnm not in ALIGNMENTS_BY_RNAME:
				ALIGNMENTS_BY_RNAME[rnm] = []
			ALIGNMENTS_BY_RNAME[rnm].append([read_pos_1, read_pos_2, ref, pos1, pos2, orientation, mapq, rdat])
		f.close()
	#
	elif INPUT_TYPE == 'pickle':
		print('reading pickle input...')
		f = open(INPUT_SAM, 'rb')
		ALIGNMENTS_BY_RNAME = pickle.load(f)
		f.close()
		sorted_keys = sorted(ALIGNMENTS_BY_RNAME.keys())

	if INPUT_TYPE != 'pickle':
		print('sorting reads...')
		sorted_keys = sorted([(sorted_keys[k][0], sorted_keys[k][1], k) for k in sorted_keys.keys()])
		sorted_keys = [n[2] for n in sorted_keys]

	#
	#
	#
	#
	#
	print('processing reads...')
	ALL_ANCHORS   = {}
	reads_skipped = {'min_readlen':0,
	                 'no_tel':0,
	                 'tel_frac':0,
	                 'only_tel':0,
	                 'unexplained_seq':0,
	                 'nontel_kmer_dens':0,
	                 'unmapped':0,
	                 'trim_filter':0}
	NONTEL_REF_SPANS_BY_CHR = {}
	#
	for job_i in range(MYJOB-1, len(sorted_keys), NJOBS):

		rnm = sorted_keys[job_i]
		if PRINT_DEBUG:
			print('read', job_i, '/', len(sorted_keys))

		abns_k = repeated_matches_trimming(sorted(ALIGNMENTS_BY_RNAME[rnm]), strategy=MATCH_TRIM_STRATEGY, print_debug=PRINT_DEBUG)

		# did we lose all of our alignments during trimming?
		if len(abns_k) == 0:
			reads_skipped['trim_filter'] += 1
			continue

		# make sure string used for kmer matching is same orientation as the alignments
		which_i = 0
		for i in range(len(abns_k)):
			if abns_k[i][2][:3] != 'tel':
				which_i = i
				break
		# assuming softclipping was used. i.e. all alignments should have same sequence... (don't pick tel though)
		if abns_k[which_i][5] == 'FWD':
			rdat = abns_k[which_i][7]
		elif abns_k[which_i][5] == 'REV':
			rdat = RC(abns_k[which_i][7])

		# read len filter
		if len(rdat) < MINIMUM_READ_LEN:
			reads_skipped['min_readlen'] += 1
			continue

		# check if we're unmapped
		refs_we_aln_to = [aln[2] for aln in abns_k]
		refs_we_aln_to = sorted(list(set(refs_we_aln_to)))

		# compute telomere kmer density
		(td_p_e0, td_p_e1) = get_telomere_kmer_density(rdat, 'p', TEL_WINDOW_SIZE)
		(td_q_e0, td_q_e1) = get_telomere_kmer_density(rdat, 'q', TEL_WINDOW_SIZE)

		#
		#	ESTIMATE TELOMERE BOUNDARIES
		#
		(p_vs_q_power, tel_regions) = get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, TEL_WINDOW_SIZE, P_VS_Q_AMP_THRESH)
		#
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
		elif which_to_choose == 'q':
			my_tel_len_p  = 0
			my_tel_len    = my_tel_len_q
			my_plot_title = ' : q-arm tel = ' + str(my_tel_len_q)
			my_tel_aln    = [[len(rdat)-my_tel_len_q, len(rdat), 'tel-q', None, None, 'REV']]

		#
		# filters
		#
		skip_me = False
		# (1) unmapped
		if skip_me == False:
			if refs_we_aln_to == ['*']:
				reads_skipped['unmapped'] += 1
				my_plot_dir = OUT_PLOT_DIR + 'filt-unmapped/'
				skip_me     = True
		# (2) too much tel
		my_tel_frac = (my_tel_len_p + my_tel_len_q) / float(len(rdat))
		if my_tel_frac > MAXIMUM_TEL_FRAC:
			reads_skipped['tel_frac'] += 1
			my_plot_dir = OUT_PLOT_DIR + 'filt-telfrac/'
			skip_me     = True
		# (3) no tel at all
		if skip_me == False:
			if my_tel_len_p == 0 and my_tel_len_q == 0:
				reads_skipped['no_tel'] += 1
				my_plot_dir = OUT_PLOT_DIR + 'filt-notel/'
				if INPUT_TYPE != 'pickle':
					for aln in abns_k:
						if aln[2] not in NONTEL_REF_SPANS_BY_CHR:
							NONTEL_REF_SPANS_BY_CHR[aln[2]] = []
						NONTEL_REF_SPANS_BY_CHR[aln[2]].append(tuple(sorted(aln[3:5])))
				skip_me = True
		# (4) too much unexplained seq
		if skip_me == False:
			aln_cov = np.zeros(len(rdat))
			for aln in abns_k:
				if aln[2][:3] != 'tel':
					aln_cov[aln[0]:aln[1]] = 1
			(p1, p2) = (my_tel_len_p, len(rdat) - my_tel_len_q)
			my_unexplained_frac = 1.0 - (np.sum(aln_cov[p1:p2]) / float(p2 - p1))
			if my_unexplained_frac > MAXIMUM_UNEXPLAINED_FRAC:
				reads_skipped['unexplained_seq'] += 1
				my_plot_dir = OUT_PLOT_DIR + 'filt-unexplained/'
				skip_me     = True
		# (5) non-tel regions look too much like tel...
		if skip_me == False:
			my_nontel_kmer_dens = abs(np.median(p_vs_q_power[my_tel_len_p:len(p_vs_q_power)-my_tel_len_q]))
			if my_nontel_kmer_dens > MAX_NONTEL_MEDIAN_KMER_DENSITY:
				reads_skipped['nontel_kmer_dens'] += 1
				my_plot_dir = OUT_PLOT_DIR + 'filt-nontelkmer/'
				skip_me     = True

		#
		#	ANCHOR TELS TO SUBTELOMERE ALIGNMENTS
		#
		if skip_me:
			pass
		else:
			#
			#	FIND ANCHOR POINT
			#
			alns_with_tel       = []
			dist_to_nearest_aln = None
			adjacent_chr        = None
			adjacent_pos        = None
			adjacent_span       = None
			my_tel_type         = None
			tel_ref_span        = None	# estimated ref coords spanned by tel
			only_tel_filter     = False
			my_plot_dir         = ''
			if my_tel_len_p > 0:
				my_tel_pos = my_tel_len_p
				for i in range(len(abns_k)):
					if abns_k[i][2][:3] == 'tel' or abns_k[i][1] <= my_tel_pos:	# tel supersedes
						continue
					if abns_k[i][0] <= my_tel_pos and abns_k[i][1] > my_tel_pos:
						alns_with_tel.append(copy.deepcopy(abns_k[i]))
						adj = my_tel_pos - abns_k[i][0]
						alns_with_tel[-1][0] += adj
						if alns_with_tel[-1][5] == 'FWD':
							alns_with_tel[-1][3] += adj
							tel_ref_span = (alns_with_tel[-1][3] - my_tel_len_p, alns_with_tel[-1][3])
						elif alns_with_tel[-1][5] == 'REV':
							alns_with_tel[-1][3] -= adj
							tel_ref_span = (alns_with_tel[-1][3], alns_with_tel[-1][3] + my_tel_len_p)
					else:
						alns_with_tel.append(copy.deepcopy(abns_k[i]))
				if len(alns_with_tel) == 0:
					reads_skipped['only_tel'] += 1
					my_plot_dir = OUT_PLOT_DIR + 'filt-onlytel/'
					only_tel_filter = True
				#
				if only_tel_filter == False:
					alns_with_tel = my_tel_aln + alns_with_tel
					if ANCHORING_STRATEGY == 'closest':
						ati = 1
					elif ANCHORING_STRATEGY == 'largest':
						ati = sorted([(alns_with_tel[n][1] - alns_with_tel[n][0], n) for n in range(len(alns_with_tel)) if alns_with_tel[n][2][:3] != 'tel'])[-1][1]
					#
					dist_to_nearest_aln = alns_with_tel[ati][0] - alns_with_tel[0][1]
					adjacent_chr        = alns_with_tel[ati][2]
					adjacent_pos        = alns_with_tel[ati][3]
					adjacent_span       = sorted([alns_with_tel[ati][4], alns_with_tel[ati][3]])
					my_tel_type         = 'p'
					if tel_ref_span == None:
						if alns_with_tel[ati][5] == 'FWD':
							tel_ref_span = (alns_with_tel[ati][3] - my_tel_len_p - dist_to_nearest_aln, alns_with_tel[ati][3] - dist_to_nearest_aln)
						elif alns_with_tel[ati][5] == 'REV':
							tel_ref_span = (alns_with_tel[ati][3] + dist_to_nearest_aln, alns_with_tel[ati][3] + my_tel_len_p + dist_to_nearest_aln)
			#
			elif my_tel_len_q > 0:
				my_tel_pos = len(rdat)-my_tel_len_q
				for i in range(len(abns_k)):
					if abns_k[i][2][:3] == 'tel' or abns_k[i][0] >= my_tel_pos:	# tel supersedes
						continue
					if abns_k[i][0] <= my_tel_pos and abns_k[i][1] > my_tel_pos:
						alns_with_tel.append(copy.deepcopy(abns_k[i]))
						adj = abns_k[i][1] - my_tel_pos
						alns_with_tel[-1][1] -= adj
						if alns_with_tel[-1][5] == 'FWD':
							alns_with_tel[-1][4] -= adj
							tel_ref_span = (alns_with_tel[-1][4], alns_with_tel[-1][4] + my_tel_len_q)
						elif alns_with_tel[-1][5] == 'REV':
							alns_with_tel[-1][4] += adj
							tel_ref_span = (alns_with_tel[-1][4] - my_tel_len_q, alns_with_tel[-1][4])
					else:
						alns_with_tel.append(copy.deepcopy(abns_k[i]))
				if len(alns_with_tel) == 0:
					reads_skipped['only_tel'] += 1
					my_plot_dir = OUT_PLOT_DIR + 'filt-onlytel/'
					only_tel_filter = True
				#
				if only_tel_filter == False:
					alns_with_tel = alns_with_tel + my_tel_aln
					if ANCHORING_STRATEGY == 'closest':
						ati = -2
					elif ANCHORING_STRATEGY == 'largest':
						ati = sorted([(alns_with_tel[n][1] - alns_with_tel[n][0], n) for n in range(len(alns_with_tel)) if alns_with_tel[n][2][:3] != 'tel'])[-1][1]
					#
					dist_to_nearest_aln = alns_with_tel[-1][0] - alns_with_tel[ati][1]
					adjacent_chr        = alns_with_tel[ati][2]
					adjacent_pos        = alns_with_tel[ati][4]
					adjacent_span       = sorted([alns_with_tel[ati][4], alns_with_tel[ati][3]])
					my_tel_type         = 'q'
					if tel_ref_span == None:
						if alns_with_tel[ati][5] == 'FWD':
							tel_ref_span = (alns_with_tel[ati][4] + dist_to_nearest_aln, alns_with_tel[ati][4] + my_tel_len_q + dist_to_nearest_aln)
						elif alns_with_tel[ati][5] == 'REV':
							tel_ref_span = (alns_with_tel[ati][4] - my_tel_len_q - dist_to_nearest_aln, alns_with_tel[ati][4] - dist_to_nearest_aln)
			#
			for aln in alns_with_tel:
				if aln[2][:3] != 'tel':
					if aln[2] not in NONTEL_REF_SPANS_BY_CHR:
						NONTEL_REF_SPANS_BY_CHR[aln[2]] = []
					NONTEL_REF_SPANS_BY_CHR[aln[2]].append(tuple(sorted(aln[3:5])))
				if PRINT_DEBUG:
					print(job_i, aln[:7])
			if PRINT_DEBUG:
				print(adjacent_chr, adjacent_pos, my_tel_len, tel_ref_span, dist_to_nearest_aln, 'ind:', ati)
			#
			if adjacent_chr not in ALL_ANCHORS:
				ALL_ANCHORS[adjacent_chr] = []
			ALL_ANCHORS[adjacent_chr].append([rnm, adjacent_pos, adjacent_span, my_tel_len, my_tel_type, tel_ref_span, copy.deepcopy(rdat), [n[:7] for n in alns_with_tel]])
			#
			if my_plot_dir == '':
				my_plot_dir = OUT_PLOT_DIR + adjacent_chr + '/'

		#	PLOTTING
		#
		if PLOT_READS:
			makedir(my_plot_dir)
			plot_title = str(job_i) + ' : ' + rnm + my_plot_title
			fig_name   = my_plot_dir + 'read_' + str(job_i) + '.png'
			dens_data  = [td_p_e0, td_p_e1, td_q_e0, td_q_e1, p_vs_q_power]
			tlen_vals  = [my_tel_len_p, my_tel_len_q]
			#
			plot_all_read_data(dens_data, tlen_vals, abns_k, TEL_WINDOW_SIZE, plot_title, fig_name)

	# sort non-tel ref spans
	#
	for k in NONTEL_REF_SPANS_BY_CHR.keys():
		NONTEL_REF_SPANS_BY_CHR[k] = sorted(NONTEL_REF_SPANS_BY_CHR[k])

	# print stats on reads that were filtered
	#
	print()
	for k in sorted(reads_skipped.keys()):
		print(k, reads_skipped[k])
	print()

	# save a pickle
	#
	f = open(OUT_PICKLE, 'wb')
	pickle.dump({'anchored-tels':ALL_ANCHORS, 'non-tel-ref-spans':NONTEL_REF_SPANS_BY_CHR}, f)
	f.close()

if __name__ == '__main__':
	main()
