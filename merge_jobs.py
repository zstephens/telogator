import argparse
import bisect
import copy
import os
import pathlib
import pickle
import re
import sys

import numpy as np

from source.tg_cluster import cluster_tel_sequences, convert_colorvec_to_kmerhits
from source.tg_kmer    import get_nonoverlapping_kmer_hits, get_telomere_kmer_frac
from source.tg_plot    import plot_kmer_hits, tel_len_violin_plot, anchor_confusion_matrix
from source.tg_util    import RC, cluster_list, LEXICO_2_IND, exists_and_is_nonzero, makedir

#
MIN_DOUBLE_ANCHOR_LEN   = 1000
MIN_DOUBLE_ANCHOR_READS = 3
#
DUMMY_TEL_MAPQ = 60
#
INCLUDE_SUBTEL_BUFF = 500	# how much into the subtel alignment should we consider when clustering reads

#
#
#
def main(raw_args=None):
	parser = argparse.ArgumentParser(description='merge_jobs.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,                  metavar='* in_dir/',      help="* Path to telogator results directory")
	parser.add_argument('-r',  type=str, required=False, default='hifi', metavar='[hifi]',         help="Type of reads (hifi, clr, ont)")
	parser.add_argument('-t',  type=str, required=False, default='p90',  metavar='[p90]',          help="Method for computing chr TL (mean/median/max/p90)")
	parser.add_argument('-cd', type=int, required=False, default=2000,   metavar='[2000]',         help="Maximum distance apart to cluster anchor positions")
	parser.add_argument('-cp', type=int, required=False, default=15,     metavar='[15]',           help="Mininum percentage of tel that must CCCTAA (15 = 15%)")
	parser.add_argument('-rc', type=int, required=False, default=2,      metavar='[2]',            help="Minimum number of reads per cluster")
	parser.add_argument('-ra', type=int, required=False, default=3,      metavar='[3]',            help="Minimum number of reads per phased allele")
	parser.add_argument('-ta', type=str, required=False, default='max',  metavar='[max]',          help="Method for computing allele TL (mean/median/max/p90)")
	parser.add_argument('-tc', type=str, required=False, default='',     metavar='treecuts.tsv',   help="Custom treecut vals during allele clustering")
	#
	parser.add_argument('-th', type=int, required=False, default=0,        metavar='[0]',            help="TelomereHunter tel_content (for comparison)")
	parser.add_argument('-gt', type=str, required=False, default='',       metavar='tlens.tsv',      help="Ground truth tel lens (for comparison)")
	parser.add_argument('-rl', type=int, required=False, default=50000,    metavar='[50000]',        help="Maximum y-axis value for readlength violin plots")
	parser.add_argument('-k',  type=str, required=False, default='',       metavar='plot_kmers.tsv', help="Telomere kmers to use for composition plotting")
	parser.add_argument('-m',  type=str, required=False, default='muscle', metavar='muscle',         help="/path/to/muscle executable")
	#
	parser.add_argument('--pbsim',                 required=False, default=False,  action='store_true', help='Simulated data from pbsim, print out confusion matrix')
	parser.add_argument('--tel-composition-plots', required=False, default=False,  action='store_true', help='Produce telomere sequence composition plots')
	parser.add_argument('--plot-denoised-tel',     required=False, default=False,  action='store_true', help='Use denoised reads for sequence composition plotting')
	parser.add_argument('--extra-tlen-plots',      required=False, default=False,  action='store_true', help='Produce extra violin plots (TL)')
	parser.add_argument('--extra-readlen-plots',   required=False, default=False,  action='store_true', help='Produce extra violin plots (readlens)')
	parser.add_argument('--telogator-pickle',      required=False, default=False,  action='store_true', help='Produce pickle which can be reprocessed by telogator.py')
	args = parser.parse_args()

	IN_DIR = args.i
	if IN_DIR[-1] != '/':
		IN_DIR += '/'
	IN_PREFIX = 'tel-data'

	READ_TYPE = args.r.lower()
	if READ_TYPE not in ['hifi', 'clr', 'ont']:
		print('Error: -m must be hifi, clr, or ont')
		exit(1)

	TL_METHOD        = args.t
	TL_METHOD_ALLELE = args.ta

	SUMMARY_SCATTER = IN_DIR + 'tlens_scatter.png'
	OUT_TSV         = IN_DIR + 'tlens_by_chr.tsv'
	MERGED_ALNS     = IN_DIR + 'merged_aln.p'
	TH_COMP_OUT     = IN_DIR + 'telomerehunter_comparison.tsv'
	CONFUSION_PLOT  = IN_DIR + 'subtel_confusion_matrix.png'

	VIOLIN_00 = IN_DIR + 'tlens_violin_00_all.png'
	VIOLIN_01 = IN_DIR + 'tlens_violin_01_chr.png'
	VIOLIN_02 = IN_DIR + 'tlens_violin_02_alt.png'
	VIOLIN_03 = IN_DIR + 'tlens_violin_03_all-pass.png'
	VIOLIN_04 = IN_DIR + 'tlens_violin_04_all-filt.png'
	VIOLIN_05 = IN_DIR + 'tlens_violin_05_chr-pass.png'
	VIOLIN_06 = IN_DIR + 'tlens_violin_06_chr-filt.png'
	VIOLIN_07 = IN_DIR + 'tlens_violin_07_alt-pass.png'
	VIOLIN_08 = IN_DIR + 'tlens_violin_08_alt-filt.png'
	#
	READLEN_VIOLIN_03 = IN_DIR + 'readlens_violin_03_all-pass.png'
	READLEN_VIOLIN_04 = IN_DIR + 'readlens_violin_04_all-filt.png'
	READLEN_VIOLIN_05 = IN_DIR + 'readlens_violin_05_chr-pass.png'
	READLEN_VIOLIN_06 = IN_DIR + 'readlens_violin_06_chr-filt.png'
	READLEN_VIOLIN_07 = IN_DIR + 'readlens_violin_07_alt-pass.png'
	READLEN_VIOLIN_08 = IN_DIR + 'readlens_violin_08_alt-filt.png'
	#
	TEL_SEQUENCES_FASTA  = IN_DIR + 'sequences-tel.fa'
	SUB_SEQUENCES_FASTA  = IN_DIR + 'sequences-sub.fa'
	READ_SEQUENCES_FASTA = IN_DIR + 'sequences-read.fa'
	#
	TELCOMP_DIR    = IN_DIR + 'phased_tel_composition/'
	DENDROGRAM_DIR = IN_DIR + 'phased_tel_dendrograms/'
	DISTMATRIX_DIR = IN_DIR + 'phased_tel_distmatrix/'
	TRAINING_DIR   = IN_DIR + 'phased_tel_dist-train/'
	CONSENSUS_DIR  = IN_DIR + 'phased_tel_msa/'
	ALLELE_OUT_FN  = IN_DIR + 'tlens_by_allele.tsv'

	ANCHOR_CLUSTER_DIST = args.cd
	MIN_READS_PER_CLUST = args.rc
	MIN_READS_PER_PHASE = args.ra
	READS_ARE_PBSIM     = args.pbsim
	TEL_HUNTER_AVG      = args.th
	GROUND_TRUTH_TLENS  = args.gt
	TEL_SEQ_PLOTS       = args.tel_composition_plots
	PLOT_DENOISED_CVECS = args.plot_denoised_tel
	EXTRA_TLEN_PLOTS    = args.extra_tlen_plots
	EXTRA_READLEN_PLOTS = args.extra_readlen_plots
	EXTRA_READLEN_YMAX  = args.rl
	TELOGATOR_PICKLE    = args.telogator_pickle
	TREECUT_TSV         = args.tc
	MIN_CANONICAL_FRAC  = args.cp/100.
	MUSCLE_EXE          = args.m

	if TEL_SEQ_PLOTS:
		makedir(TELCOMP_DIR)
		makedir(DENDROGRAM_DIR)
		makedir(DISTMATRIX_DIR)
		makedir(TRAINING_DIR)
		makedir(CONSENSUS_DIR)

	KMER_FILE   = args.k
	KMER_LIST   = []
	KMER_COLORS = []
	rev_kmers   = []
	rev_colors  = []
	if KMER_FILE == '':
		print('using default telomere kmers.')
		sim_path = pathlib.Path(__file__).resolve().parent
		kmer_fn  = str(sim_path) + '/resources/plot_kmers.tsv'
		f = open(kmer_fn,'r')
		for line in f:
			if line[0] != '#':
				splt = line.strip().split('\t')
				KMER_LIST.append(splt[1])
				KMER_COLORS.append(splt[2])
				rev_kmers.append(RC(splt[1]))
				rev_colors.append(splt[2])
		f.close()
	else:
		fn_suffix = KMER_FILE.split('/')[-1]
		print('using user-specified kmer list:', fn_suffix)
		if exists_and_is_nonzero(KMER_FILE):
			f = open(KMER_FILE,'r')
			for line in f:
				if line[0] != '#':
					splt = line.strip().split('\t')
					KMER_LIST.append(splt[1])
					KMER_COLORS.append(splt[2])
					rev_kmers.append(RC(splt[1]))
					rev_colors.append(splt[2])
			f.close()
		else:
			print('Error: kmer list not found')
			exit(1)
	sorted_kmer_dat  = sorted(list(set([(len(KMER_LIST[n]),KMER_LIST[n],KMER_COLORS[n]) for n in range(len(KMER_LIST))])), reverse=True)	# sort by length
	KMER_LIST        = [n[1] for n in sorted_kmer_dat]
	KMER_COLORS      = [n[2] for n in sorted_kmer_dat]
	KMER_ISSUBSTRING = []
	for i in range(len(KMER_LIST)):
		KMER_ISSUBSTRING.append([j for j in range(len(KMER_LIST)) if (j != i and KMER_LIST[i] in KMER_LIST[j])])
		#print(i, KMER_LIST[i], KMER_COLORS[i], KMER_ISSUBSTRING[-1])
	KMER_LIST_REV = [RC(n) for n in KMER_LIST]

	#
	# read in ground truth values for plotting (e.g. for simulated data)
	#
	gt_tlen = {}
	if len(GROUND_TRUTH_TLENS):
		f = open(GROUND_TRUTH_TLENS, 'r')
		for line in f:
			splt = line.strip().split('\t')
			gt_tlen[splt[0].replace('_','')] = int(splt[2])
		f.close()

	#
	# cuts for specific chr in hand-picked places
	#
	custom_treecut_vals = {}
	if len(TREECUT_TSV):
		f = open(TREECUT_TSV, 'r')
		for line in f:
			if len(line.strip()):
				splt = line.strip().split('\t')
				if ':' in splt[0]:
					splt2  = splt[0].split(':')
					tc_chr = splt2[0]
					tc_pos = int(splt2[1])
				else:
					tc_chr = splt[0].replace('_','')
					tc_pos = None
				custom_treecut_vals[(tc_chr, tc_pos)] = float(splt[1])
		f.close()

	#
	# look for jobs
	#
	listing = [n for n in os.listdir(IN_DIR) if (n[:len(IN_PREFIX)] == IN_PREFIX and '_job-' in n and n[-2:] == '.p')]
	if len(listing) == 0:
		print('Error: no pickles found in input dir matching prefix:', IN_PREFIX)
		exit(1)

	sorted_fns = {}
	for fn in listing:
		splt   = fn[:-2].split('_')[-1].split('-')
		my_key = int(splt[2])
		if my_key not in sorted_fns:
			sorted_fns[my_key] = []
		sorted_fns[my_key].append((int(splt[1]), fn))
	if len(sorted_fns) != 1:
		print('Error: jobs with different run totals found with same prefix')
		exit(1)
	sorted_fns = [n[1] for n in sorted(sorted_fns[my_key])]

	#
	COMBINED_ANCHORS  = {}
	COMBINED_REFSPANS = {}
	#
	for fn in sorted_fns:
		f = open(IN_DIR + fn, 'rb')
		my_pickle = pickle.load(f)
		f.close()
		#
		my_anchors = my_pickle['anchored-tels']
		for k in my_anchors.keys():
			if k not in COMBINED_ANCHORS:
				COMBINED_ANCHORS[k] = []
			COMBINED_ANCHORS[k].extend(copy.deepcopy(my_anchors[k]))
		#
		my_refspans = my_pickle['non-tel-ref-spans']
		for k in my_refspans.keys():
			if k not in COMBINED_REFSPANS:
				COMBINED_REFSPANS[k] = []
			COMBINED_REFSPANS[k].extend(copy.deepcopy(my_refspans[k]))

	# need to resort after job merge
	for k in sorted(COMBINED_REFSPANS.keys()):
		COMBINED_REFSPANS[k] = sorted(COMBINED_REFSPANS[k])

	# ANCHORED-TEL LIST FORMAT:
	#
	#  0    1             2              3           4            5             6     7
	# [rnm, adjacent_pos, adjacent_span, my_tel_len, my_tel_type, tel_ref_span, rdat, alns_with_tel]

	#
	# filter out double-anchored telomeres
	#
	print('filtering out double-anchored tel regions...')
	#
	del_keys = []
	for k in sorted(COMBINED_ANCHORS.keys()):
		len_before_filtering = len(COMBINED_ANCHORS[k])
		del_list = []
		if k in COMBINED_REFSPANS:
			for i in range(len(COMBINED_ANCHORS[k])):
				tel_ref_span   = COMBINED_ANCHORS[k][i][5]
				relevant_spans = []
				####print(COMBINED_ANCHORS[k][i][:6], tel_ref_span)
				bi  = bisect.bisect(COMBINED_REFSPANS[k], tel_ref_span)
				bi2 = min([bi, len(COMBINED_REFSPANS[k]) - 1])
				while True:
					if bi2 < 0 or tel_ref_span[0] - COMBINED_REFSPANS[k][bi2][1] > 10*MIN_DOUBLE_ANCHOR_LEN:
						break
					relevant_spans.append(COMBINED_REFSPANS[k][bi2])
					bi2 -= 1
				bi2 = min([bi, len(COMBINED_REFSPANS[k]) - 1])
				while True:
					if bi2 >= len(COMBINED_REFSPANS[k]) or COMBINED_REFSPANS[k][bi2][0] - tel_ref_span[1] > 10*MIN_DOUBLE_ANCHOR_LEN:
						break
					relevant_spans.append(COMBINED_REFSPANS[k][bi2])
					bi2 += 1
				relevant_spans = sorted(relevant_spans)
				n_spanning_reads = 0
				for span in relevant_spans:
					if tel_ref_span[0] - span[0] >= MIN_DOUBLE_ANCHOR_LEN and span[1] - tel_ref_span[1] >= MIN_DOUBLE_ANCHOR_LEN:
						n_spanning_reads += 1
				if n_spanning_reads >= MIN_DOUBLE_ANCHOR_READS:
					del_list.append(i)
		del_list = sorted(del_list, reverse=True)
		for di in del_list:
			del COMBINED_ANCHORS[k][di]
		del del_list
		#print(k, len_before_filtering, '-->', len(COMBINED_ANCHORS[k]))
		if len(COMBINED_ANCHORS[k]) == 0:
			del_keys.append(k)
	for k in del_keys:
		del COMBINED_ANCHORS[k]

	#
	ANCHORED_TEL_ALL  = [{}, {}, {}]
	ANCHORED_TEL_FILT = [{}, {}, {}, {}, {}, {}]
	READLEN_FILT      = [{}, {}, {}, {}, {}, {}]
	CHR_OR_ALT        = {}
	CONF_DAT          = {}
	sorted_ref_keys   = []
	#
	for k in sorted(COMBINED_ANCHORS.keys()):
		if k[:3] == 'chr':
			my_chr = k.replace('_','')
			my_ind = 1
		elif k[:3] == 'alt':
			my_chr = ''.join(k.split('_')[:-1])
			my_chr = my_chr.replace('alt', 'chr')
			my_ind = 2
		else:
			print('skipping ref:', k)
			continue
		#
		sorted_ref_keys.append((my_ind, LEXICO_2_IND[my_chr[:-1]], my_chr, k))
		#
		if my_chr not in ANCHORED_TEL_ALL[0]:
			ANCHORED_TEL_ALL[0][my_chr] = []
			ANCHORED_TEL_ALL[1][my_chr] = []
			ANCHORED_TEL_ALL[2][my_chr] = []
		for i in range(len(COMBINED_ANCHORS[k])):
			ANCHORED_TEL_ALL[0][my_chr].append(COMBINED_ANCHORS[k][i][3])
			ANCHORED_TEL_ALL[my_ind][my_chr].append(COMBINED_ANCHORS[k][i][3])
			#
			rnm = COMBINED_ANCHORS[k][i][0]
			if my_ind == 1:
				CHR_OR_ALT[rnm] = 'chr'
			elif my_ind == 2:
				CHR_OR_ALT[rnm] = 'alt'
			#
			if READS_ARE_PBSIM:
				conf_key = (rnm.split('-')[1].replace('_',''), my_chr)
				if conf_key not in CONF_DAT:
					CONF_DAT[conf_key] = 0
				CONF_DAT[conf_key] += 1
	sorted_ref_keys = sorted(sorted_ref_keys)

	# reconstruct reduced aln set for plotting/debugging:
	#
	# [readpos_start, readpos_end, ref, pos_start, pos_end, orientation, mapq, read_sequence]
	#
	COMBINED_ALN_DAT = {}

	#
	#
	#
	CHR_TEL_DAT = []
	#
	if TEL_SEQ_PLOTS:
		f_telfasta = open(TEL_SEQUENCES_FASTA, 'w')
		f_subfasta = open(SUB_SEQUENCES_FASTA, 'w')
		f_entire   = open(READ_SEQUENCES_FASTA, 'w')
		ALLELE_CLUST_DAT = []
	#
	unexplained_telseq_dict = {}
	allele_count_by_chr     = {}
	#
	len_by_read_tel = []
	len_by_read_sub = []
	#
	comp_data = []
	clust_num = 0
	for ki in range(len(sorted_ref_keys)):
		k = sorted_ref_keys[ki][3]
		#
		if k[:3] == 'chr':
			my_chr = k.replace('_','')
		elif k[:3] == 'alt':
			my_chr = ''.join(k.split('_')[:-1])
			my_chr = my_chr.replace('alt', 'chr')
		print(k, my_chr, len(COMBINED_ANCHORS[k]))
		#
		sort_list = sorted([(COMBINED_ANCHORS[k][n][1], n) for n in range(len(COMBINED_ANCHORS[k]))])
		clusters  = cluster_list(sort_list, ANCHOR_CLUSTER_DIST, which_val=0)
		#
		for cl in clusters:
			pos_list = [n[0] for n in cl]
			ind_list = [n[1] for n in cl]

			#
			# remove reads with tel regions that do not contain a large enough fraction of canonical repeats
			#
			del_list = []
			for di in range(len(ind_list)):
				my_tlen = COMBINED_ANCHORS[k][ind_list[di]][3]
				my_type = COMBINED_ANCHORS[k][ind_list[di]][4]
				my_rdat = COMBINED_ANCHORS[k][ind_list[di]][6]
				if my_chr[-1] == 'p':
					if my_type == 'p':
						my_telseq = my_rdat[:my_tlen]
					elif my_type == 'q':
						my_telseq = RC(my_rdat[-my_tlen:])
				elif my_chr[-1] == 'q':
					if my_type == 'p':
						my_telseq = RC(my_rdat[:my_tlen])
					elif my_type == 'q':
						my_telseq = my_rdat[-my_tlen:]
				my_canonical_frac = get_telomere_kmer_frac(my_telseq, ['CCCTAA', 'TTAGGG'], mode=READ_TYPE)
				if my_canonical_frac < MIN_CANONICAL_FRAC:
					del_list.append(di)
			del_list = sorted(del_list, reverse=True)
			for di in del_list:
				del pos_list[di]
				del ind_list[di]
			del del_list

			#
			# determine which violin plots we will be contributing our data to
			#
			sorted_by_tlen = sorted([(COMBINED_ANCHORS[k][i][3], len(COMBINED_ANCHORS[k][i][6]), COMBINED_ANCHORS[k][i][0]) for i in ind_list])
			my_tlens       = [n_sbtl[0] for n_sbtl in sorted_by_tlen]
			my_rlens       = [n_sbtl[1] for n_sbtl in sorted_by_tlen]
			my_rnames      = [n_sbtl[2] for n_sbtl in sorted_by_tlen]
			for i in range(len(my_rnames)):
				my_rname = my_rnames[i]
				if CHR_OR_ALT[my_rname] == 'chr' and len(cl) >= MIN_READS_PER_CLUST:
					my_inds = [0, 2]
				elif CHR_OR_ALT[my_rname] == 'chr' and len(cl) < MIN_READS_PER_CLUST:
					my_inds = [1, 3]
				elif CHR_OR_ALT[my_rname] == 'alt' and len(cl) >= MIN_READS_PER_CLUST:
					my_inds = [0, 4]
				elif CHR_OR_ALT[my_rname] == 'alt' and len(cl) < MIN_READS_PER_CLUST:
					my_inds = [1, 5]
				for my_ind in my_inds:
					if my_chr not in ANCHORED_TEL_FILT[my_ind]:
						ANCHORED_TEL_FILT[my_ind][my_chr] = []
						READLEN_FILT[my_ind][my_chr]      = []
					ANCHORED_TEL_FILT[my_ind][my_chr].append(my_tlens[i])
					READLEN_FILT[my_ind][my_chr].append(my_rlens[i])

			#
			# if the cluster does not enough reads lets skip it
			#
			if len(ind_list) < MIN_READS_PER_CLUST or len(ind_list) <= 0:
				continue
			my_pos = int(np.median(pos_list))
			clust_num += 1

			#
			# output combined data for easier reprocessing by telogator.py if desired (e.g. for plotting individual reads / kmer densities)
			#
			for i in ind_list:
				my_rnm  = COMBINED_ANCHORS[k][i][0]
				my_rdat = COMBINED_ANCHORS[k][i][6]
				my_alns = COMBINED_ANCHORS[k][i][7]
				# revert RC using same logic as telogator
				which_i = 0
				for i in range(len(my_alns)):
					if my_alns[i][2][:3] != 'tel':
						which_i = i
						break
				if my_alns[which_i][5] == 'REV':
					my_rdat = RC(my_rdat)
				#
				if my_rnm not in COMBINED_ALN_DAT:
					COMBINED_ALN_DAT[my_rnm] = []
				for aln in my_alns:
					COMBINED_ALN_DAT[my_rnm].append([copy.deepcopy(n) for n in aln] + [DUMMY_TEL_MAPQ]*(aln[2][:3] == 'tel') + [my_rdat])

			#print('---', len(pos_list), int(np.mean(pos_list)), int(np.std(pos_list)))
			#for i in ind_list:
			#	print('-------', COMBINED_ANCHORS[k][i][:5])

			consensus_tl = None
			if TL_METHOD == 'mean':
				consensus_tl = np.mean(my_tlens)
			elif TL_METHOD == 'median':
				consensus_tl = np.median(my_tlens)
			elif TL_METHOD == 'max':
				consensus_tl = np.max(my_tlens)
			elif TL_METHOD[0] == 'p':
				my_percentile = int(TL_METHOD[1:])
				consensus_tl = np.percentile(my_tlens, my_percentile)

			CHR_TEL_DAT.append([sorted_ref_keys[ki][3],
			                    str(my_pos),
			                    str(int(consensus_tl)),
			                    ','.join([str(n) for n in my_tlens]),
			                    ','.join([str(n) for n in my_rlens])])
			#
			comp_data.append([sorted_ref_keys[ki][3], my_pos, [n for n in my_tlens]])

			########################################################
			#
			# plot kmer composition of telomeres in this cluster
			#
			########################################################
			if TEL_SEQ_PLOTS:
				kmer_hit_dat = []
				rlens_out    = []
				for i in ind_list:
					my_rnm  = COMBINED_ANCHORS[k][i][0]
					my_tlen = COMBINED_ANCHORS[k][i][3]
					my_type = COMBINED_ANCHORS[k][i][4]
					my_rdat = COMBINED_ANCHORS[k][i][6]
					my_alns = COMBINED_ANCHORS[k][i][7]
					#
					rlens_out.append(len(my_rdat))
					#
					my_mapq    = None
					my_dbta    = None	# distance between telomere and anchor
					anchor_mai = None
					my_anchor_ref_coords = sorted(COMBINED_ANCHORS[k][i][2])
					for mai in range(len(my_alns)):
						my_aln = my_alns[mai]
						if my_aln[3]!= None and my_aln[4] != None and sorted([my_aln[3], my_aln[4]]) == my_anchor_ref_coords:
							anchor_mai = mai
					if anchor_mai != None:
						my_mapq = my_alns[anchor_mai][6]
						if my_alns[0][2][:3] == 'tel':		# tel is first alignment
							my_dbta = my_alns[anchor_mai][0] - my_alns[0][1]
						elif my_alns[-1][2][:3] == 'tel':	# tel is last alignment
							my_dbta = my_alns[-1][0] - my_alns[anchor_mai][1]
						if my_dbta == None or my_dbta < 0:
							print('Error: we messed up trying to get your telomere-anchor dist:', k, i, my_dbta)
					#
					# adjusted telomere boundary
					#
					atb = my_tlen + my_dbta + INCLUDE_SUBTEL_BUFF
					atb = min(atb, len(my_rdat))					# in case something crazy happens
					#
					if my_chr[-1] == 'p':
						kmers_to_use = KMER_LIST
						if my_type == 'p':
							my_telseq = my_rdat[:atb]
							my_subseq = my_rdat[atb:]
						elif my_type == 'q':
							my_telseq = RC(my_rdat[-atb:])
							my_subseq = RC(my_rdat[:-atb])
					elif my_chr[-1] == 'q':
						kmers_to_use = KMER_LIST_REV
						if my_type == 'p':
							my_telseq = RC(my_rdat[:atb])
							my_subseq = RC(my_rdat[atb:])
						elif my_type == 'q':
							my_telseq = my_rdat[-atb:]
							my_subseq = my_rdat[:-atb]
					len_by_read_tel.append(len(my_telseq))
					len_by_read_sub.append(len(my_subseq))
					#
					out_readname = 'cluster-' + str(clust_num) + '_ref-' + my_chr + '_tel-' + my_type + '_' + my_rnm
					f_telfasta.write('>' + out_readname + '\n')
					f_telfasta.write(my_telseq + '\n')
					out_readname = 'cluster-' + str(clust_num) + '_ref-' + my_chr + '_sub-' + my_type + '_' + my_rnm
					f_subfasta.write('>' + out_readname + '\n')
					f_subfasta.write(my_subseq + '\n')
					out_readname = my_rnm
					f_entire.write('>' + out_readname + '\n')
					f_entire.write(my_rdat + '\n')
					#
					# get kmer hits
					#
					# kmer_hit_dat[-1][0][kmer_list_i] = hits in current read for kmer kmer_list_i
					#
					kmer_hit_dat.append([get_nonoverlapping_kmer_hits(my_telseq, kmers_to_use, KMER_ISSUBSTRING), atb, my_dbta, my_type, my_rnm, my_mapq])
					#
					# what are the unexplained sequences?
					#
					coord_hit_all = np.zeros(atb)
					for kmer_list_i in range(len(KMER_LIST)):
						for kmer_span in kmer_hit_dat[-1][0][kmer_list_i]:
							for j in range(kmer_span[0],kmer_span[1]):
								coord_hit_all[j] = 1
					unexplained_regions = []
					for j in range(atb):
						if coord_hit_all[j] == 0:
							unexplained_regions.append([j,j+1])
					for j in range(len(unexplained_regions)-1,0,-1):
						if unexplained_regions[j-1][1] == unexplained_regions[j][0]:
							unexplained_regions[j-1][1] = unexplained_regions[j][1]
							del unexplained_regions[j]
					for j in range(len(unexplained_regions)):
						ur = unexplained_regions[j]
						us = my_telseq[ur[0]:ur[1]]
						if us not in unexplained_telseq_dict:
							unexplained_telseq_dict[us] = 0
						unexplained_telseq_dict[us] += 1
				#
				# plotting!
				#
				if k[:3] == 'alt':
					plotname_chr = 'alt-' + my_chr
				else:
					plotname_chr = my_chr
				#
				if True or plotname_chr == 'chr3p':
					zfcn = str(clust_num).zfill(2)
					dendrogram_fn  = DENDROGRAM_DIR + 'cluster-' + zfcn + '_' + plotname_chr + '.png'
					distmatrix_fn  = DISTMATRIX_DIR + 'cluster-' + zfcn + '_' + plotname_chr + '.npy'
					telcompplot_fn = TELCOMP_DIR    + 'cluster-' + zfcn + '_' + plotname_chr + '.png'
					telcompcons_fn = TELCOMP_DIR    + 'msa-'     + zfcn + '_' + plotname_chr + '.png'
					traindata_pref = TRAINING_DIR   + 'cluster-' + zfcn + '_' + plotname_chr
					consensus_fn   = CONSENSUS_DIR  + 'cluster-' + zfcn + '_' + plotname_chr + '.fa'
					#
					my_tc = None
					if (my_chr,None) in custom_treecut_vals:
						my_tc = custom_treecut_vals[(my_chr,None)]
					elif (my_chr,my_pos) in custom_treecut_vals:
						my_tc = custom_treecut_vals[(my_chr,my_pos)]
					#
					read_clust_dat = cluster_tel_sequences(kmer_hit_dat, KMER_LIST, KMER_COLORS, my_chr, my_pos, dist_in=distmatrix_fn, fig_name=dendrogram_fn, tree_cut=my_tc, msa_dir=IN_DIR, train_prefix=traindata_pref, save_msa=consensus_fn, muscle_exe=MUSCLE_EXE)
					#
					for allele_i in range(len(read_clust_dat[0])):
						#
						allele_tvrlen        = read_clust_dat[7][allele_i]
						allele_consensus_out = ''
						if allele_tvrlen > 0:
							if my_chr[-1] == 'p':	# p will be reversed so its in subtel --> tvr --> tel orientation
								allele_consensus_out = read_clust_dat[4][allele_i][-allele_tvrlen:][::-1]
							elif my_chr[-1] == 'q':
								allele_consensus_out = read_clust_dat[4][allele_i][:allele_tvrlen]
						#
						# kmer_hit_dat[n][1]   = tlen + all the extra subtel buffers
						# read_clust_dat[3][n] = the length of the subtel region present before tvr/tel region
						#
						# the difference of these two will be the actual size of the (tvr + tel) region in the read
						#
						allele_readcount = len(read_clust_dat[0][allele_i])
						allele_tlen_mapq = sorted([(kmer_hit_dat[n][1] - read_clust_dat[3][n], rlens_out[n], kmer_hit_dat[n][5]) for n in read_clust_dat[0][allele_i]])
						allele_tlens     = [n[0]-len(allele_consensus_out) for n in allele_tlen_mapq]	# subtracting tvr so that "actual" TL is output. values can be negative
						allele_tlen_str  = ','.join([str(n) for n in allele_tlens])
						rlen_str         = ','.join([str(n[1]) for n in allele_tlen_mapq])
						mapq_str         = ','.join([str(n[2]) for n in allele_tlen_mapq])
						#
						consensus_tl_allele = None
						if TL_METHOD_ALLELE == 'mean':
							consensus_tl_allele = np.mean(allele_tlens)
						elif TL_METHOD_ALLELE == 'median':
							consensus_tl_allele = np.median(allele_tlens)
						elif TL_METHOD_ALLELE == 'max':
							consensus_tl_allele = np.max(allele_tlens)
						elif TL_METHOD_ALLELE[0] == 'p':
							my_percentile = int(TL_METHOD_ALLELE[1:])
							consensus_tl_allele = np.percentile(allele_tlens, my_percentile)
						#
						if allele_readcount >= MIN_READS_PER_PHASE:
							if my_chr not in allele_count_by_chr:
								allele_count_by_chr[my_chr] = 0
							ALLELE_CLUST_DAT.append([my_chr,
								                    str(my_pos),
								                    str(allele_count_by_chr[my_chr]),
								                    str(int(consensus_tl_allele)),
								                    allele_tlen_str,
								                    rlen_str,
								                    mapq_str,
								                    str(len(allele_consensus_out)),
								                    allele_consensus_out])
							allele_count_by_chr[my_chr] += 1
					#
					# adjust kmer_hit_dat based on the filters and etc that were applied during clustering
					#
					my_consensus_vecs = read_clust_dat[4]
					my_color_vectors  = read_clust_dat[5]
					my_end_err_lens   = read_clust_dat[6]
					my_tvr_tel_bounds = read_clust_dat[7]
					redrawn_kmerhits  = convert_colorvec_to_kmerhits(my_color_vectors, KMER_LIST, KMER_COLORS)
					redrawn_consensus = convert_colorvec_to_kmerhits(my_consensus_vecs, KMER_LIST, KMER_COLORS)
					if PLOT_DENOISED_CVECS:
						for rdki in range(len(redrawn_kmerhits)):
							kmer_hit_dat[rdki][0]  = redrawn_kmerhits[rdki]	# replace kmer hit tuples for plotting
							kmer_hit_dat[rdki][1] -= my_end_err_lens[rdki]	# subtract size of artifacts at end of reads
					#
					consensus_kmer_hit_dat = []
					consensus_clust_dat    = [[],[],[],[0]]	# fake data so that plot_kmer_hits doesn't complain
					consensus_tvr_tel_pos  = []
					for rdki in range(len(redrawn_consensus)):
						cons_readcount = len(read_clust_dat[0][rdki])
						cons_readname  = 'consensus-' + str(rdki) + ' [' + str(cons_readcount)
						cons_tvrlen    = my_tvr_tel_bounds[rdki]
						if cons_readcount == 1:
							cons_readname += ' read]'
						else:
							cons_readname += ' reads]'
						if cons_readcount >= MIN_READS_PER_PHASE:
							consensus_kmer_hit_dat.append([redrawn_consensus[rdki], len(my_consensus_vecs[rdki]), 0, 'FWD', cons_readname, DUMMY_TEL_MAPQ])
							consensus_clust_dat[0].append([rdki])
							consensus_clust_dat[1].append([DUMMY_TEL_MAPQ])
							consensus_clust_dat[2].append([0])
							consensus_tvr_tel_pos.append(cons_tvrlen)
					#
					# plotting!
					#
					plot_kmer_hits(kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompplot_fn, xlim=[-1000,15000], clust_dat=read_clust_dat)
					if len(consensus_clust_dat[0]):
						plot_kmer_hits(consensus_kmer_hit_dat, KMER_COLORS, my_chr, my_pos, telcompcons_fn, xlim=[-1000,15000], clust_dat=consensus_clust_dat, draw_boundaries=consensus_tvr_tel_pos)
		print()
	#
	# write chr-level output (old telogator)
	#
	f_chr = open(OUT_TSV, 'w')
	f_chr.write('#chr' + '\t' + 'position' + '\t' + 'TL_' + TL_METHOD + '\t' + 'read_TLs' + '\t' + 'read_lengths' + '\n')
	for i in range(len(CHR_TEL_DAT)):
		f_chr.write('\t'.join(CHR_TEL_DAT[i]) + '\n')
	f_chr.close()
	#
	# write allele-level output
	#
	if TEL_SEQ_PLOTS:
		f_entire.close()
		f_subfasta.close()
		f_telfasta.close()
		#
		f_allele = open(ALLELE_OUT_FN, 'w')
		f_allele.write('#chr' + '\t' + 'position' + '\t')
		f_allele.write('allele_num' + '\t' + 'allele_TL_' + TL_METHOD_ALLELE + '\t' + 'read_allele_TLs' + '\t' + 'read_lengths' + '\t' + 'read_MAPQ' + '\t' + 'TVR_length' + '\t' + 'TVR_consensus' + '\n')
		for i in range(len(ALLELE_CLUST_DAT)):
			f_allele.write('\t'.join(ALLELE_CLUST_DAT[i]) + '\n')
		f_allele.close()
		#
		os.system('cat ' + TRAINING_DIR + '*_tvrtel.tsv > ' + TRAINING_DIR + 'all-tvrtel.tsv')

	#
	if TELOGATOR_PICKLE:
		f = open(MERGED_ALNS, 'wb')
		pickle.dump(COMBINED_ALN_DAT, f)
		f.close()

	my_pm = True
	if len(gt_tlen):
		my_pm = False

	readlen_plot_params = {'p_color':'gray',
	                       'q_color':'gray',
	                       'y_label':'<-- q     read length     p -->',
	                       'p_ymax':EXTRA_READLEN_YMAX,
	                       'q_ymax':EXTRA_READLEN_YMAX,
	                       'y_step':10000}

	if EXTRA_TLEN_PLOTS:
		tel_len_violin_plot(ANCHORED_TEL_ALL[0], VIOLIN_00, plot_means=my_pm, ground_truth_dict=gt_tlen)
		tel_len_violin_plot(ANCHORED_TEL_ALL[1], VIOLIN_01, plot_means=my_pm, ground_truth_dict=gt_tlen)
		tel_len_violin_plot(ANCHORED_TEL_ALL[2], VIOLIN_02, plot_means=my_pm, ground_truth_dict=gt_tlen)
	#
	tel_len_violin_plot(ANCHORED_TEL_FILT[0], VIOLIN_03, plot_means=my_pm, ground_truth_dict=gt_tlen)
	if EXTRA_READLEN_PLOTS:
		tel_len_violin_plot(READLEN_FILT[0], READLEN_VIOLIN_03, custom_plot_params=readlen_plot_params)
	if EXTRA_TLEN_PLOTS:
		tel_len_violin_plot(ANCHORED_TEL_FILT[1], VIOLIN_04, plot_means=my_pm, ground_truth_dict=gt_tlen)
		if EXTRA_READLEN_PLOTS:
			tel_len_violin_plot(READLEN_FILT[1], READLEN_VIOLIN_04, custom_plot_params=readlen_plot_params)
	#
	tel_len_violin_plot(ANCHORED_TEL_FILT[2], VIOLIN_05, plot_means=my_pm, ground_truth_dict=gt_tlen)
	if EXTRA_READLEN_PLOTS:
		tel_len_violin_plot(READLEN_FILT[2], READLEN_VIOLIN_05, custom_plot_params=readlen_plot_params)
	if EXTRA_TLEN_PLOTS:
		tel_len_violin_plot(ANCHORED_TEL_FILT[3], VIOLIN_06, plot_means=my_pm, ground_truth_dict=gt_tlen)
		if EXTRA_READLEN_PLOTS:
			tel_len_violin_plot(READLEN_FILT[3], READLEN_VIOLIN_06, custom_plot_params=readlen_plot_params)
	#
	tel_len_violin_plot(ANCHORED_TEL_FILT[4], VIOLIN_07, plot_means=my_pm, ground_truth_dict=gt_tlen)
	if EXTRA_READLEN_PLOTS:
		tel_len_violin_plot(READLEN_FILT[4], READLEN_VIOLIN_07, custom_plot_params=readlen_plot_params)
	if EXTRA_TLEN_PLOTS:
		tel_len_violin_plot(ANCHORED_TEL_FILT[5], VIOLIN_08, plot_means=my_pm, ground_truth_dict=gt_tlen)
		if EXTRA_READLEN_PLOTS:
			tel_len_violin_plot(READLEN_FILT[5], READLEN_VIOLIN_08, custom_plot_params=readlen_plot_params)

	if READS_ARE_PBSIM:
		anchor_confusion_matrix(CONF_DAT, CONFUSION_PLOT)

	####if TEL_SEQ_PLOTS:
	####	skus = sorted([(unexplained_telseq_dict[k],k) for k in unexplained_telseq_dict.keys()], reverse=True)
	####	print('unexplained tel regions:')
	####	for n in skus:
	####		if n[0] >= 2:
	####			print(n[0], n[1])
	####	print('')

	#
	if TEL_HUNTER_AVG > 0 and len(comp_data):
		print('stats for telomerehunter comparison...')
		f = open(TH_COMP_OUT, 'w')
		f.write('tel_content' + '\t' + str(TEL_HUNTER_AVG) + '\n')
		comp_name  = ['mean', 'median', 'max', 'p90']
		#
		print('- all')
		comp_print = [[], [], [], []]
		for n in comp_data:
			comp_print[0].append(np.mean(n[2]))
			comp_print[1].append(np.median(n[2]))
			comp_print[2].append(max(n[2]))
			comp_print[3].append(np.percentile(n[2],90))
		for i in range(len(comp_name)):
			print('---', comp_name[i], ':', np.mean(comp_print[i]))
			f.write('all' + '\t' + comp_name[i] + '\t' + str(str(int(np.mean(comp_print[i])))) + '\n')
		#
		print('- chr')
		comp_print = [[], [], [], []]
		for n in comp_data:
			if n[0][:3] == 'chr':
				comp_print[0].append(np.mean(n[2]))
				comp_print[1].append(np.median(n[2]))
				comp_print[2].append(max(n[2]))
				comp_print[3].append(np.percentile(n[2],90))
		for i in range(len(comp_name)):
			print('---', comp_name[i], ':', np.mean(comp_print[i]))
			f.write('chr' + '\t' + comp_name[i] + '\t' + str(str(int(np.mean(comp_print[i])))) + '\n')
		f.close()


if __name__ == '__main__':
	main()
