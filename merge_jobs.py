import argparse
import bisect
import copy
import os
import pathlib
import pickle
import re
import sys

import numpy as np

from source.tg_cluster import cluster_tel_sequences
from source.tg_kmer    import get_nonoverlapping_kmer_hits
from source.tg_plot    import plot_kmer_hits, tel_len_violin_plot, anchor_confusion_matrix
from source.tg_util    import RC, cluster_list, LEXICO_2_IND, exists_and_is_nonzero

#
MIN_DOUBLE_ANCHOR_LEN   = 1000
MIN_DOUBLE_ANCHOR_READS = 3
#
DUMMY_TEL_MAPQ = 60
#

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='merge_jobs.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,                  metavar='* in_dir/',      help="* Path to telogator results directory")
	parser.add_argument('-t',  type=str, required=False, default='p90',  metavar='[p90]',          help="Method for computing TL (mean/median/max/p90)")
	parser.add_argument('-cd', type=int, required=False, default=2000,   metavar='[2000]',         help="Maximum distance apart to cluster anchor positions")
	parser.add_argument('-cr', type=int, required=False, default=2,      metavar='[2]',            help="Minimum number of reads per cluster")
	parser.add_argument('-th', type=int, required=False, default=0,      metavar='[0]',            help="TelomereHunter tel_content (for comparison)")
	parser.add_argument('-gt', type=str, required=False, default='',     metavar='tlens.tsv',      help="Ground truth tel lens (for comparison)")
	parser.add_argument('-rl', type=int, required=False, default=50000,  metavar='[50000]',        help="Maximum y-axis value for readlength violin plots")
	parser.add_argument('-k',  type=str, required=False, default='',     metavar='plot_kmers.tsv', help="Telomere kmers to use for composition plotting")
	#
	parser.add_argument('--pbsim',                 required=False, default=False,  action='store_true', help='Simulated data from pbsim, print out confusion matrix')
	parser.add_argument('--tel-composition-plots', required=False, default=False,  action='store_true', help='Produce telomere sequence composition plots')
	parser.add_argument('--extra-tlen-plots',      required=False, default=False,  action='store_true', help='Produce extra violin plots (TL)')
	parser.add_argument('--extra-readlen-plots',   required=False, default=False,  action='store_true', help='Produce extra violin plots (readlens)')
	args = parser.parse_args()

	IN_DIR = args.i
	if IN_DIR[-1] != '/':
		IN_DIR += '/'
	IN_PREFIX = 'tel-data'

	TL_METHOD = args.t

	SUMMARY_SCATTER = IN_DIR + 'tel_lens_scatter.png'
	CONFUSION_PLOT  = IN_DIR + 'subtel_confusion_matrix.png'
	OUT_TSV         = IN_DIR + 'results.tsv'
	MERGED_ALNS     = IN_DIR + 'merged_aln.p'
	TH_COMP_OUT     = IN_DIR + 'telomerehunter_comparison.tsv'

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
	TEL_SEQUENCES_FASTA  = IN_DIR + 'tel-sequences.fa'
	SUB_SEQUENCES_FASTA  = IN_DIR + 'sub-sequences.fa'
	READ_SEQUENCES_FASTA = IN_DIR + 'read-sequences.fa'

	ANCHOR_CLUSTER_DIST = args.cd
	MIN_READS_PER_CLUST = args.cr
	READS_ARE_PBSIM     = args.pbsim
	TEL_HUNTER_AVG      = args.th
	GROUND_TRUTH_TLENS  = args.gt
	TEL_SEQ_PLOTS       = args.tel_composition_plots
	EXTRA_TLEN_PLOTS    = args.extra_tlen_plots
	EXTRA_READLEN_PLOTS = args.extra_readlen_plots
	EXTRA_READLEN_YMAX  = args.rl

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
			splt = line.strip().split('\t')
			KMER_LIST.append(splt[1])
			KMER_COLORS.append(splt[2])
			rev_kmers.append(RC(splt[1]))
			rev_colors.append(splt[2])
		f.close()
		KMER_LIST.extend([n for n in rev_kmers])
		KMER_COLORS.extend([n for n in rev_colors])
	else:
		fn_suffix = KMER_FILE.split('/')[-1]
		print('using user-specified kmer list:', fn_suffix)
		if exists_and_is_nonzero(KMER_FILE):
			f = open(KMER_FILE,'r')
			for line in f:
				splt = line.strip().split('\t')
				KMER_LIST.append(splt[1])
				KMER_COLORS.append(splt[2])
				rev_kmers.append(RC(splt[1]))
				rev_colors.append(splt[2])
			f.close()
			KMER_LIST.extend([n for n in rev_kmers])
			KMER_COLORS.extend([n for n in rev_colors])
		else:
			print('Error: kmer list not found')
			exit(1)
	sorted_kmer_dat  = sorted([(len(KMER_LIST[n]),KMER_LIST[n],KMER_COLORS[n]) for n in range(len(KMER_LIST))], reverse=True)	# sort by length
	KMER_LIST        = [n[1] for n in sorted_kmer_dat]
	KMER_COLORS      = [n[2] for n in sorted_kmer_dat]
	KMER_ISSUBSTRING = []
	for i in range(len(KMER_LIST)):
		KMER_ISSUBSTRING.append([j for j in range(len(KMER_LIST)) if (j != i and KMER_LIST[i] in KMER_LIST[j])])
	####for i in range(len(KMER_LIST)):
	####	print(i, KMER_LIST[i], KMER_COLORS[i], KMER_ISSUBSTRING[i])
	####exit(1)

	gt_tlen = {}
	if len(GROUND_TRUTH_TLENS):
		f = open(GROUND_TRUTH_TLENS, 'r')
		for line in f:
			splt = line.strip().split('\t')
			gt_tlen[splt[0].replace('_','')] = int(splt[2])
		f.close()

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
			fail = False
			for n in my_anchors[k]:
				if n[3] > 15000:
					fail = True
			if fail:
				continue
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

	# filter out double-anchored telomeres
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

	####for k in sorted(ANCHORED_TEL_ALL[0].keys()):
	####	print(k)
	####	print(ANCHORED_TEL_ALL[0][k])
	####	print()

	sorted_ref_keys = sorted(sorted_ref_keys)

	# reconstruct reduced aln set for plotting/debugging:
	#
	# [readpos_start, readpos_end, ref, pos_start, pos_end, orientation, mapq, read_sequence]
	#
	COMBINED_ALN_DAT = {}

	#
	f_out = open(OUT_TSV, 'w')
	f_out.write('#subtel' + '\t' + 'position' + '\t' + 'tel_len_' + TL_METHOD + '\t' + 'tel_lens' + '\t' + 'read_lens' + '\n')
	#
	if TEL_SEQ_PLOTS:
		f_telfasta = open(TEL_SEQUENCES_FASTA, 'w')
		f_subfasta = open(SUB_SEQUENCES_FASTA, 'w')
		f_entire   = open(READ_SEQUENCES_FASTA, 'w')
	#
	unexplained_telseq_dict = {}
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
		for cl in clusters:
			#
			clust_num += 1
			pos_list   = [n[0] for n in cl]
			ind_list   = [n[1] for n in cl]
			my_pos     = int(np.mean(pos_list))
			#
			#my_rnames = [COMBINED_ANCHORS[k][i][0] for i in ind_list]
			#my_tlens  = [COMBINED_ANCHORS[k][i][3] for i in ind_list]
			#my_rlens  = [len(COMBINED_ANCHORS[k][i][6]) for i in ind_list]
			#
			sorted_by_tlen = sorted([(COMBINED_ANCHORS[k][i][3], len(COMBINED_ANCHORS[k][i][6]), COMBINED_ANCHORS[k][i][0]) for i in ind_list])
			my_tlens  = [n_sbtl[0] for n_sbtl in sorted_by_tlen]
			my_rlens  = [n_sbtl[1] for n_sbtl in sorted_by_tlen]
			my_rnames = [n_sbtl[2] for n_sbtl in sorted_by_tlen]

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

			if len(cl) < MIN_READS_PER_CLUST:
				continue

			#
			# output combined data for easier reprocessing if desired
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

			print('---', len(pos_list), int(np.mean(pos_list)), int(np.std(pos_list)))
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

			f_out.write(sorted_ref_keys[ki][3] + '\t' + str(my_pos) + '\t' + str(int(consensus_tl)) + '\t' + ','.join([str(n) for n in my_tlens]) + '\t' + ','.join([str(n) for n in my_rlens]) + '\n')
			#
			comp_data.append([sorted_ref_keys[ki][3], my_pos, [n for n in my_tlens]])

			#
			# plot kmer composition of telomeres in this cluster
			#
			if TEL_SEQ_PLOTS:
				kmer_hit_dat    = []
				for i in ind_list:
					my_rnm  = COMBINED_ANCHORS[k][i][0]
					my_tlen = COMBINED_ANCHORS[k][i][3]
					my_type = COMBINED_ANCHORS[k][i][4]
					my_rdat = COMBINED_ANCHORS[k][i][6]
					my_alns = COMBINED_ANCHORS[k][i][7]
					#print('OUTPUT TEL FASTA:', my_chr, COMBINED_ANCHORS[k][i][:6])
					if my_chr[-1] == 'p':
						if my_type == 'p':
							my_telseq = my_rdat[:my_tlen]
							my_subseq = my_rdat[my_tlen:]
						elif my_type == 'q':
							my_telseq = RC(my_rdat[-my_tlen:])
							my_subseq = RC(my_rdat[:-my_tlen])
					elif my_chr[-1] == 'q':
						if my_type == 'p':
							my_telseq = RC(my_rdat[:my_tlen])
							my_subseq = RC(my_rdat[my_tlen:])
						elif my_type == 'q':
							my_telseq = my_rdat[-my_tlen:]
							my_subseq = my_rdat[:-my_tlen]
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
					kmer_hit_dat.append([get_nonoverlapping_kmer_hits(my_telseq, KMER_LIST, KMER_ISSUBSTRING), my_tlen, my_type, my_rnm])
					#
					# what are the unexplained sequences?
					#
					coord_hit_all = np.zeros(my_tlen)
					for kmer_list_i in range(len(KMER_LIST)):
						for kmer_span in kmer_hit_dat[-1][0][kmer_list_i]:
							for j in range(kmer_span[0],kmer_span[1]):
								coord_hit_all[j] = 1
					unexplained_regions = []
					for j in range(my_tlen):
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
				#####
				####dendrogram_fn = IN_DIR + 'dendrogram_cluster-' + str(clust_num) + '_ref-' + plotname_chr + '.png'
				####read_clust_dat = cluster_tel_sequences(kmer_hit_dat, KMER_COLORS, my_chr, fig_name=dendrogram_fn)
				#####
				####plot_fn = IN_DIR + 'tel-composition-plot_cluster-' + str(clust_num) + '_ref-' + plotname_chr + '.png'
				####plot_kmer_hits(kmer_hit_dat, KMER_COLORS, my_chr, plot_fn, clust_dat=read_clust_dat)
		print()
	#
	if TEL_SEQ_PLOTS:
		f_entire.close()
		f_subfasta.close()
		f_telfasta.close()
		#
		####import matplotlib.pyplot as mpl
		####mpl.figure(0)
		####mpl.scatter(len_by_read_sub, len_by_read_tel)
		####mpl.show()
		####exit(1)
	f_out.close()

	#
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

	if TEL_SEQ_PLOTS:
		skus = sorted([(unexplained_telseq_dict[k],k) for k in unexplained_telseq_dict.keys()], reverse=True)
		print('unexplained tel regions:')
		for n in skus:
			if n[0] >= 2:
				print(n[0], n[1])
		print('')

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
