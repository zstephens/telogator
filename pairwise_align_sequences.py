import os
import sys

import numpy as np
import matplotlib.pyplot as mpl
import matplotlib.colors as colors

from source.tg_cluster import cluster_consensus_tel
from source.tg_reader import TG_Reader
from source.tg_util import makedir

def main(raw_args=None):
	#
	in_fasta = sys.argv[1]
	out_dir  = sys.argv[2]
	num_proc = int(sys.argv[3])
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
	out_dist = out_dir + 'dist-matrix.npy'
	out_fig  = out_dir + 'dendrogram.png'
	#
	cluster_consensus_tel(all_tvrs,
		                  dist_in=out_dist,
		                  fig_name=out_fig,
		                  samp_labels=my_labels,
		                  tree_cut=0.20,
		                  alignment_processes=num_proc)

if __name__ == '__main__':
	main()
