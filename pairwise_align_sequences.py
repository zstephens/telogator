import argparse
import os
import pathlib
import sys

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors

from source.tg_cluster import cluster_consensus_tel, convert_colorvec_to_kmerhits, MIN_MSD
from source.tg_plot    import plot_kmer_hits, tel_len_bar_plot
from source.tg_reader  import TG_Reader
from source.tg_util    import exists_and_is_nonzero, makedir, RC

def main(raw_args=None):
	#
	parser = argparse.ArgumentParser(description='pairwise_align_sequences.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str,   required=True,  metavar='* in.fasta', help="* Fasta of TVRs to align")
	parser.add_argument('-o',  type=str,   required=True,  metavar='* out_dir/', help="Output directory")
	parser.add_argument('-m',  type=str,   required=False, metavar='[ms/ds]',    help="Which pairwise2 call to use (ms or ds)",   default='ms')
	parser.add_argument('-t',  type=float, required=False, metavar='0.50',       help="Dendrogram tree cut position",             default=0.50)
	parser.add_argument('-d',  type=int,   required=False, metavar='12',         help="Dendrogram figure size (height)",          default=12)
	parser.add_argument('-c',  type=int,   required=False, metavar='2',          help="Minimum #TVR alleles per cluster to plot", default=3)
	#
	parser.add_argument('-p',             type=int, required=False, metavar='8',       default=8,     help="Number of alignment processes per job")
	parser.add_argument('--job', nargs=2, type=int, required=False, metavar=('i','N'), default=(1,1), help='Jobs IDs (i/N) for running in parallel (e.g. 1 5, 9 9). i=0 to merge.')
	args = parser.parse_args()
	#
	in_fasta  = args.i
	out_dir   = args.o
	aln_mode  = args.m
	num_proc  = args.p
	tree_cut  = args.t
	dend_size = args.d
	min_allele_plot = args.c
	(my_job, total_jobs) = args.job
	#
	if out_dir[-1] != '/':
		out_dir += '/'
	makedir(out_dir)
	#
	if aln_mode not in ['ms','ds']:
		print('Error: -m must be ms or ds')
		exit(1)
	#
	#
	#
	KMER_LIST   = []
	KMER_COLORS = []
	sim_path = pathlib.Path(__file__).resolve().parent
	kmer_fn  = str(sim_path) + '/resources/plot_kmers.tsv'
	f = open(kmer_fn,'r')
	for line in f:
		if line[0] != '#':
			splt = line.strip().split('\t')
			KMER_LIST.append(splt[1])
			KMER_COLORS.append(splt[2])
	f.close()
	sorted_kmer_dat  = sorted(list(set([(len(KMER_LIST[n]),KMER_LIST[n],KMER_COLORS[n]) for n in range(len(KMER_LIST))])), reverse=True)	# sort by length
	KMER_LIST        = [n[1] for n in sorted_kmer_dat]
	KMER_COLORS      = [n[2] for n in sorted_kmer_dat]
	KMER_ISSUBSTRING = []
	for i in range(len(KMER_LIST)):
		KMER_ISSUBSTRING.append([j for j in range(len(KMER_LIST)) if (j != i and KMER_LIST[i] in KMER_LIST[j])])
	KMER_LIST_REV = [RC(n) for n in KMER_LIST]
	#
	#
	#
	my_labels = []
	all_tvrs  = []
	my_reader = TG_Reader(in_fasta, verbose=False)
	while True:
		read_dat = my_reader.get_next_read()
		if not read_dat[0]:
			break
		my_labels.append(read_dat[0])
		all_tvrs.append(read_dat[1])
	my_reader.close()
	#
	####all_tvrs  = all_tvrs[:20]	# truncation for testing
	####my_labels = my_labels[:20]	# truncation for testing
	#
	out_dist = out_dir + 'dist-matrix.npy'
	out_fig  = out_dir + 'dendrogram.png'
	#
	if my_job == 0:
		listing = sorted([n for n in os.listdir(out_dir) if ('_job' in n and n[-4:] == '.npy')])
		if len(listing) < total_jobs:
			print('Error: missing a job!')
			exit(1)
		for i in range(len(listing)):
			if i == 0:
				dist_matrix = np.load(out_dir + listing[i], allow_pickle=True)
			else:
				next_matrix = np.load(out_dir + listing[i], allow_pickle=True)
				dist_matrix = dist_matrix + next_matrix
		dist_norm    = max(np.max(dist_matrix), MIN_MSD)
		dist_matrix /= dist_norm
		np.save(out_dist, dist_matrix)
	#
	if exists_and_is_nonzero(out_dist):
		conf_dist = []
		conf_hits = {}
		conf_indv = {}
		allele_is_multimapped = {n:False for n in my_labels}
		conf_max  = 0.08
		dist_matrix = np.load(out_dist, allow_pickle=True)
		n_tvrs = dist_matrix.shape[0]
		for i in range(n_tvrs):
			for j in range(i+1,n_tvrs):
				splt_i = my_labels[i].split('_')
				splt_j = my_labels[j].split('_')
				if splt_i[0] == splt_j[0] and splt_i[1] != splt_j[1] and dist_matrix[i,j] <= conf_max:
					#print(i, j, my_labels[i], my_labels[j], dist_matrix[i,j])
					conf_dist.append(dist_matrix[i,j])
					my_key = tuple(sorted([splt_i[1], splt_j[1]]))
					if my_key not in conf_hits:
						conf_hits[my_key] = 0
					conf_hits[my_key] += 1
					if splt_i[1] not in conf_indv:
						conf_indv[splt_i[1]] = 0
					conf_indv[splt_i[1]] += 1
					if splt_j[1] not in conf_indv:
						conf_indv[splt_j[1]] = 0
					conf_indv[splt_j[1]] += 1
					allele_is_multimapped[my_labels[i]] = True
					allele_is_multimapped[my_labels[j]] = True
		sk = sorted([[conf_hits[k],k] for k in conf_hits.keys()], reverse=True)
		total_multimapped_alleles = sum(conf_hits.values())
		print('TOTAL SIMILAR TVRS:', total_multimapped_alleles, '/', n_tvrs, '({0:0.2f}%)'.format(100.*(float(total_multimapped_alleles)/n_tvrs)))
		for n in sk:
			if n[0] >= 2:
				print(n[0], '({0:0.2f}%)'.format(100.*(float(n[0])/total_multimapped_alleles)), n[1])
		sk = sorted([[conf_indv[k],k] for k in conf_indv.keys()], reverse=True)
		barplot_dat = {}
		for n in sk:
			numer = 0
			denom = 0
			for k in allele_is_multimapped.keys():
				splt = k.split('_')
				if splt[1] == n[1]:
					denom += 1
					if allele_is_multimapped[k] == True:
						numer += 1
			if denom >= 10:
				print(numer, '/', denom, '({0:0.2f}%)'.format(100.*(float(numer)/denom)), n[1])
				barplot_dat[n[1]] = 100.*(float(numer)/denom)
		#
		allelecount_dat = {}
		multimapped_dat = {}
		for k in allele_is_multimapped.keys():
			splt = k.split('_')
			if splt[1] not in allelecount_dat:
				allelecount_dat[splt[1]] = 0
				multimapped_dat[splt[1]] = 0
			allelecount_dat[splt[1]] += 1
			if allele_is_multimapped[k] == True:
				multimapped_dat[splt[1]] += 1
		#
		multimap_allele_by_chr_by_samp = {}
		for n in my_labels:
			splt = n.split('_')
			if splt[1] not in multimap_allele_by_chr_by_samp:
				multimap_allele_by_chr_by_samp[splt[1]] = {}
			if splt[0] not in multimap_allele_by_chr_by_samp[splt[1]]:
				multimap_allele_by_chr_by_samp[splt[1]][splt[0]] = [0,0]
			if allele_is_multimapped[n] == True:
				multimap_allele_by_chr_by_samp[splt[1]][splt[0]][0] += 1
			multimap_allele_by_chr_by_samp[splt[1]][splt[0]][1] += 1
		multimap_frac_by_chr = {}
		for k in multimap_allele_by_chr_by_samp.keys():
			multimap_frac_by_chr[k] = [[],[]]
			for s in multimap_allele_by_chr_by_samp[k].keys():
				multimap_frac_by_chr[k][0].append(multimap_allele_by_chr_by_samp[k][s][0])
				multimap_frac_by_chr[k][1].append(multimap_allele_by_chr_by_samp[k][s][1])
			multimap_frac_by_chr[k] = [np.mean(multimap_frac_by_chr[k][0]), np.mean(multimap_frac_by_chr[k][1])]
		mean_allele_numer = {k:multimap_frac_by_chr[k][0] for k in multimap_frac_by_chr.keys()}
		mean_allele_denom = {k:multimap_frac_by_chr[k][1] for k in multimap_frac_by_chr.keys()}
		#
		####mpl.figure(0,figsize=(10,6))
		####mpl.hist(conf_dist,  bins=50, range=[0,conf_max])
		####mpl.xlim([0,conf_max])
		####mpl.tight_layout()
		####mpl.show()
		#
		barplot_fn = out_dir + 'multimap_bar.png'
		barplot_param = {'p_ymax':100, 'q_ymax':100, 'y_step':20, 'p_color':'gray', 'q_color':'gray', 'y_label':'<-- q   multimapped alleles   p -->', 'ytick_suffix':'%'}
		tel_len_bar_plot(barplot_dat, barplot_fn, custom_plot_params=barplot_param)
		#
		barplot_fn = out_dir + 'allelecount_bar.png'
		barplot_param = {'p_ymax':300, 'q_ymax':300, 'y_step':50, 'p_color':'gray', 'q_color':'gray', 'y_label':'<-- q   allele count   p -->', 'hatch_data':multimapped_dat}
		tel_len_bar_plot(allelecount_dat, barplot_fn, custom_plot_params=barplot_param)
		#
		barplot_fn = out_dir + 'allelecount_mean_bar.png'
		barplot_param = {'p_ymax':6, 'q_ymax':6, 'y_step':1, 'p_color':'gray', 'q_color':'gray', 'y_label':'<-- q   #alleles / sample   p -->', 'hatch_data':mean_allele_numer}
		tel_len_bar_plot(mean_allele_denom, barplot_fn, custom_plot_params=barplot_param)
	exit(1)
	#
	clustdat = cluster_consensus_tel(all_tvrs,
		                             dist_in=out_dist,
		                             fig_name=out_fig,
		                             samp_labels=my_labels,
		                             aln_mode=aln_mode,
		                             tree_cut=tree_cut,
		                             alignment_processes=num_proc,
		                             job=(my_job,total_jobs),
		                             dendrogram_height=dend_size)
	#
	# plots for paper (trio)
	#
	paper_plots = [['hg002_chr12q_1', 'hg003_chr12q_0', 'hg002_chr12q_0', 'hg004_chr12q_0'],
	               ['hg002_chr2p_0', 'hg003_chr2p_1', 'hg002_chr2p_1', 'hg004_chr2p_0']]
	paper_title = ['chr12q', 'chr2p']
	paper_xlims = [(0,2000), (0,2000)]
	paper_xstep = [500, 500]
	paper_remap = {'hg002_chr12q_1':'child allele 1',
	               'hg003_chr12q_0':'father',
	               'hg002_chr12q_0':'child allele 2',
	               'hg004_chr12q_0':'mother',
	               'hg002_chr2p_0':'child allele 1',
	               'hg003_chr2p_1':'father',
	               'hg002_chr2p_1':'child allele 2',
	               'hg004_chr2p_0':'mother'}
	#
	# plots for paper (pangenome)
	#
	paper_plots = [['hg002_chr12q_1', 'hg003_chr12q_0', 'hg002_chr12q_0', 'hg004_chr12q_0']]
	paper_title = ['']
	paper_xlims = [(0,3000)]
	#
	for i in range(len(clustdat)):
		if len(clustdat[i]) >= min_allele_plot:
			print('cluster', i, len(clustdat[i]))
			clust_cvecs = [all_tvrs[ci] for ci in clustdat[i]]
			clust_label = [my_labels[ci] for ci in clustdat[i]]
			clust_kmerhits = convert_colorvec_to_kmerhits(clust_cvecs, KMER_LIST, KMER_COLORS)
			#
			# fake values needed for plotting
			#
			my_tlens = [len(all_tvrs[ci]) for ci in clustdat[i]]
			clust_kmerdat = [[clust_kmerhits[n], my_tlens[n], 0, 'q', clust_label[n], 60] for n in range(len(clust_kmerhits))]
			#
			my_chr  = 'None-q'
			my_pos  = ''
			plot_fn = out_dir + 'cluster_' + str(i+1).zfill(3) + '.png'
			plot_kmer_hits(clust_kmerdat, KMER_COLORS, my_chr, my_pos, plot_fn)
			#
			paper_kmerdat = []
			paper_fignum  = None
			for j in range(len(paper_plots)):
				found_any = False
				for n in clust_label:
					if n in paper_plots[j]:
						found_any    = True
						paper_fignum = j
						break
				if found_any:
					break
			if paper_fignum != None:
				paper_cd = [(paper_plots[paper_fignum].index(clust_label[j]), clustdat[i][j]) for j in range(len(clustdat[i])) if clust_label[j] in paper_plots[paper_fignum]]
				paper_cd = [n[1] for n in sorted(paper_cd)]
				paper_cvecs = [all_tvrs[n] for n in paper_cd]
				paper_label = []
				for n in paper_cd:
					if my_labels[n] in paper_remap:
						paper_label.append(paper_remap[my_labels[n]])
					else:
						paper_label.append(my_labels[n])
				paper_tlens = [len(all_tvrs[n]) for n in paper_cd]
				paper_kmerhits = convert_colorvec_to_kmerhits(paper_cvecs, KMER_LIST, KMER_COLORS)
				paper_kmerdat  = [[paper_kmerhits[n], paper_tlens[n], 0, 'q', paper_label[n], 60] for n in range(len(paper_kmerhits))]
				#
				my_chr  = 'None-q'
				my_pos  = ''
				paper_fn     = out_dir + 'paperfig_' + str(paper_fignum+1).zfill(3) + '.png'
				paper_params = {'custom_title':'',
				                'xlim':paper_xlims[paper_fignum],
				                'xstep':paper_xstep[paper_fignum],
				                'number_label_rows':False,
				                'custom_xlabel':'TVR region length (bp)'}
				plot_kmer_hits(paper_kmerdat, KMER_COLORS, my_chr, my_pos, paper_fn, plot_params=paper_params)

if __name__ == '__main__':
	main()
