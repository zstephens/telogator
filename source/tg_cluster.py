import os

import numpy as np
import matplotlib.pyplot as mpl

from Bio import pairwise2

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform

from source.tg_reader import TG_Reader

MUSCLE_EXE = '/Users/zach/opt/miniconda2/bin/muscle'

# we're going to pretend each kmer color is an amino acid, so alignment tools behave themselves
AMINO = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

# how many blocks of gaps to consider for plotting adj
MAX_GAP_BLOCK = 3

def write_scoring_matrix(fn):
	f = open(fn, 'w')
	f.write('   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  J  Z  X  *' + '\n')
	f.write('A  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('R -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('N -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('D -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('C -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Q -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('E -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('G -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('H -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('I -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('L -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('K -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('M -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('F -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('P -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('S -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('T -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('W -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('Y -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -4 -1 -4' + '\n')
	f.write('V -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -4 -1 -4' + '\n')
	f.write('B -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -4 -1 -4' + '\n')
	f.write('J -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -4 -1 -4' + '\n')
	f.write('Z -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  5 -1 -4' + '\n')
	f.write('X -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -4' + '\n')
	f.write('* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1' + '\n')
	f.close()

def get_muscle_msa(input_sequences, working_dir=''):
	# write sequences to a temp fasta
	temp_fasta = working_dir + 'clust_sequences.fa'
	f = open(temp_fasta, 'w')
	for i in range(len(input_sequences)):
		f.write('>seq'+str(i+1).zfill(5) + '\n')
		f.write(input_sequences[i] + '\n')
	f.close()
	# run muscle
	aln_fasta   = working_dir + 'clust_aln.fa'
	muscle_log  = working_dir + 'muscle.log'
	matrix      = working_dir + 'scoring_matrix.txt'
	write_scoring_matrix(matrix)
	score_param = '-seqtype protein -gapopen -4.0 -gapextend -1.0 -center 0.0 -matrix ' + matrix
	cmd = MUSCLE_EXE + ' -in ' + temp_fasta + ' -out ' + aln_fasta + ' ' + score_param + ' > ' + muscle_log + ' 2>&1'
	os.system(cmd)
	# get results
	out_seq   = []
	my_reader = TG_Reader(aln_fasta, verbose=False)
	while True:
		read_dat = my_reader.get_next_read()
		if not read_dat[0]:
			break
		my_readnum = int(read_dat[0][3:])
		out_seq.append([my_readnum, read_dat[1]])
	my_reader.close()
	out_seq = [n[1] for n in sorted(out_seq)]
	# cleanup
	os.system('rm ' + temp_fasta)
	os.system('rm ' + aln_fasta)
	os.system('rm ' + muscle_log)
	os.system('rm ' + matrix)
	# get consensus
	consensus_seq = []
	for i in range(len(out_seq[0])):	# we're gonna hope that muscle worked as expected and all seq are same len
		char_count = {}
		for j in range(len(out_seq)):
			if out_seq[j][i] == '-':
				continue
			if out_seq[j][i] not in char_count:
				char_count[out_seq[j][i]] = 0
			char_count[out_seq[j][i]] += 1
		sk = sorted([(char_count[k],k) for k in char_count.keys()], reverse=True)
		consensus_seq.append(sk[0][1])
	consensus_seq = ''.join(consensus_seq)
	#
	return [out_seq, consensus_seq]

def get_adj_from_gaps(s):
	in_gap     = False
	gap_count  = 0
	gap_blocks = []
	for i in range(len(s)):
		if s[i] == '-':
			if in_gap == False:
				gap_start = i
			in_gap = True
		else:
			if in_gap:
				gap_blocks.append((gap_start,i))
				gap_count += 1
			in_gap = False
		if gap_count >= MAX_GAP_BLOCK:
			break
	out_adj = 0
	if len(gap_blocks):
		gap_so_far  = [gap_blocks[0][1]-gap_blocks[0][0]]
		cost_so_far = [gap_blocks[0][0]]
		for i in range(1,len(gap_blocks)):
			gap_so_far.append(gap_so_far[-1] + gap_blocks[i][1]-gap_blocks[i][0])
			cost_so_far.append(cost_so_far[-1] + gap_blocks[i][0]-gap_blocks[i-1][1])
		best_gap = sorted([(gap_so_far[i]-cost_so_far[i],i) for i in range(len(gap_blocks))], reverse=True)[0]
		if best_gap[0] > 0:
			out_adj = gap_blocks[best_gap[1]][1]
	return out_adj

#
#	kmer_dat[i] = [[kmer1_hits, kmer2_hits, ...], tlen, read-orientation, readname, anchor_mapq]
#
def cluster_tel_sequences(kmer_dat, kmer_colors, my_chr, fig_name=None, cluster_region=2000, tree_cut=0.35, canonical_letter='C'):
	n_reads   = len(kmer_dat)
	scolors   = sorted(list(set(kmer_colors)))
	col_2_sc  = {n:scolors.index(n) for n in scolors}
	n_colors  = len(scolors)
	pq        = my_chr[-1]
	#
	# create color vector
	#
	my_col_single = []
	for i in range(n_reads):
		[my_kmer_hits, my_tlen, my_orr, my_rname, my_mapq] = kmer_dat[i]
		my_col_single.append(np.zeros(my_tlen, dtype='>i4'))
		for ki in range(len(my_kmer_hits)):
			if len(my_kmer_hits[ki]):
				for kmer_span in my_kmer_hits[ki]:
					xp = [kmer_span[0], kmer_span[1]]
					my_col_single[-1][xp[0]:xp[1]] = col_2_sc[kmer_colors[ki]]+1
		if pq == 'p':
			my_col_single[-1] = ''.join([AMINO[n] for n in my_col_single[-1][-cluster_region:]])
			pw2_gap = (False, True)
		elif pq == 'q':
			my_col_single[-1] = ''.join([AMINO[n] for n in my_col_single[-1][:cluster_region]])
			pw2_gap = (True, False)
	#
	# trivial case
	#
	if n_reads == 1:
		return [ [[0]], [[kmer_dat[0][4]]], [[0]], [my_col_single[0]]]
	#elif n_reads == 2:
	#	return [ [[0,1]], , [[0,0]] ]
	#
	# scoring matrix
	#
	letters = AMINO[:n_colors+1]
	scoring_matrix = {}
	for i in range(len(letters)):
		for j in range(len(letters)):
			if i == j:
				scoring_matrix[(letters[i],letters[j])] = 5
			else:
				scoring_matrix[(letters[i],letters[j])] = -4
	scoring_matrix[(canonical_letter, canonical_letter)] = 0
	#
	dist_matrix = np.zeros((n_reads,n_reads))
	for i in range(n_reads):
		for j in range(i+1,n_reads):
			alignments = pairwise2.align.globalms(my_col_single[i], my_col_single[j], 5, -4, -4, -4, penalize_end_gaps=pw2_gap)
			#alignments = pairwise2.align.globalds(my_col_single[i], my_col_single[j], scoring_matrix, -4, -4, penalize_end_gaps=pw2_gap)
			my_score   = int(alignments[0].score)
			aln_symbol = pairwise2.format_alignment(*alignments[0]).split('\n')[1]
			my_dist    = float(aln_symbol.count('.') + aln_symbol.count(' '))/cluster_region
			dist_matrix[i,j] = my_dist
			dist_matrix[j,i] = my_dist
			#print(pairwise2.format_alignment(*alignments[0]))
			#print(i, j, 'pairwise-distance: {0:0.3f}'.format(my_dist), 'aln-score:', my_score)
	#
	dist_array = squareform(dist_matrix)
	Zread      = linkage(dist_array, method='average')
	#
	if fig_name != None:
		fig = mpl.figure(3, figsize=(8,8))
		dendrogram(Zread)
		mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
		mpl.title(my_chr)
		mpl.savefig(fig_name)
		mpl.close(fig)
	#
	assignments = fcluster(Zread, tree_cut, 'distance').tolist()
	by_class = {}
	for i in range(len(assignments)):
		if assignments[i] not in by_class:
			by_class[assignments[i]] = []
		by_class[assignments[i]].append(i)
	out_clust = sorted([(len(by_class[k]), sorted(by_class[k])) for k in by_class.keys()], reverse=True)
	out_clust = [n[1] for n in out_clust]
	#
	out_adj       = []
	out_consensus = []
	for i in range(len(out_clust)):
		if len(out_clust[i]) == 1:
			out_adj.append([0])
			out_consensus.append(my_col_single[out_clust[i][0]])
		else:
			clust_seq                = [my_col_single[n] for n in out_clust[i]]
			[msa_seq, consensus_seq] = get_muscle_msa(clust_seq)
			out_adj.append([])
			for j in range(len(msa_seq)):
				if pq == 'q':
					out_adj[-1].append(get_adj_from_gaps(msa_seq[j]))
				elif pq == 'p':
					out_adj[-1].append(get_adj_from_gaps(msa_seq[j][::-1]))
			out_consensus.append(consensus_seq)
	#
	out_mapq = []
	for n in out_clust:
		out_mapq.append([kmer_dat[m][4] for m in n])
	#
	print(out_clust)
	print(out_mapq)
	print(out_adj)
	#print(out_consensus)
	#
	return [ out_clust, out_mapq, out_adj, out_consensus ]

#
#
#
def cluster_consensus_tel(sequences, fig_name=None):
	n_seq = len(sequences)
	dist_matrix = np.zeros((n_seq,n_seq))
	pw2_gap = (False, False)
	for i in range(n_seq):
		for j in range(i+1,n_seq):
			alignments = pairwise2.align.globalms(sequences[i], sequences[j], 5, -4, -4, -4, penalize_end_gaps=pw2_gap)
			#alignments = pairwise2.align.globalds(sequences[i], sequences[j], scoring_matrix, -4, -4, penalize_end_gaps=False)
			my_score   = int(alignments[0].score)
			aln_symbol = pairwise2.format_alignment(*alignments[0]).split('\n')[1]
			my_dist    = float(aln_symbol.count('.') + aln_symbol.count(' '))
			dist_matrix[i,j] = my_dist
			dist_matrix[j,i] = my_dist
			#print(pairwise2.format_alignment(*alignments[0]))
			print(i, j, 'pairwise-distance: {0:0.3f}'.format(my_dist), 'aln-score:', my_score)
	#
	dist_array = squareform(dist_matrix)
	Zread      = linkage(dist_array, method='average')
	#
	if fig_name != None:
		fig = mpl.figure(3, figsize=(8,8))
		dendrogram(Zread)
		#mpl.axhline(y=[tree_cut], linestyle='dashed', color='black', alpha=0.7)
		mpl.title(my_chr)
		mpl.savefig(fig_name)
		mpl.close(fig)
