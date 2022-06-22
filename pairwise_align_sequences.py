import os
import sys

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors

from source.tg_cluster import cluster_consensus_tel, MIN_MSD
from source.tg_reader import TG_Reader
from source.tg_util import makedir

def main(raw_args=None):
	#
	in_fasta = sys.argv[1]
	out_dir  = sys.argv[2]
	num_proc = int(sys.argv[3])
	(my_job, total_jobs) = (int(sys.argv[4]), int(sys.argv[5]))
	#
	if out_dir[-1] != '/':
		out_dir += '/'
	makedir(out_dir)
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
	cluster_consensus_tel(all_tvrs,
		                  dist_in=out_dist,
		                  fig_name=out_fig,
		                  samp_labels=my_labels,
		                  tree_cut=0.20,
		                  alignment_processes=num_proc,
		                  job=(my_job,total_jobs))

if __name__ == '__main__':
	main()
