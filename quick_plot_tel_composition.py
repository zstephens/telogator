import argparse
import bisect
import copy
import os
import pathlib
import pickle
import re
import sys

import numpy as np

from source.tg_plot import plot_kmer_hits
from source.tg_reader import TG_Reader
from source.tg_util import RC, cluster_list, LEXICO_2_IND

DUMMY_MAPQ = 60	# for fasta input

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='quick_plot_tel_composition.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i',  type=str, required=True,                  metavar='* input.sam',    help="* Input SAM or FASTA")
	parser.add_argument('-o',  type=str, required=True,                  metavar='* output.png',   help="* Output plot")
	parser.add_argument('-k',  type=str, required=False, default='',     metavar='plot_kmers.tsv', help="Telomere kmers to use for composition plotting")
	args = parser.parse_args()

	IN_SAM   = args.i
	OUT_PLOT = args.o

	KMER_FILE   = args.k
	KMER_LIST   = []
	KMER_COLORS = []
	if KMER_FILE == '':
		print('using default telomere kmers.')
		sim_path = pathlib.Path(__file__).resolve().parent
		kmer_fn  = str(sim_path) + '/resources/plot_kmers.tsv'
		f = open(kmer_fn,'r')
		rev_kmers  = []
		rev_colors = []
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
	
	read_dat = []
	if IN_SAM[-4:].lower() == '.sam':
		f_sam = open(IN_SAM, 'r')
		fr    = f_sam.read()
		f_sam.close()
		for line in fr.split('\n'):
			splt = line.strip().split('\t')
			read_dat.append([splt[0], len(splt[9]), splt[9], int(splt[4])])
	else:
		my_reader = TG_Reader(IN_SAM)
		while True:
			(my_name, my_rdat, my_qdat) = my_reader.get_next_read()
			if not my_name:
				break
			read_dat.append([my_name, len(my_rdat), my_rdat, DUMMY_MAPQ])
		my_reader.close()
	#
	my_chr       = 'placeholder-title q'
	my_type      = 'q'
	kmer_hit_dat = []
	line_num     = 0
	#
	for rd in read_dat:
		(my_rnm, my_tlen, my_rdat, my_mapq) = rd
		#
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
		#
		####if line_num%2 == 1:
		####	my_subseq = RC(my_subseq)
		####	my_telseq = RC(my_telseq)
		line_num += 1
		#
		kmer_hit_dat.append([[], my_tlen, my_type, my_rnm, my_mapq])	# kmer_hit_dat[-1][0][ki] = hits in current read for kmer ki
		coord_hit_dict = []	# this is sloppy but easy
		coord_hit_all  = np.zeros(my_tlen)
		for ki in range(len(KMER_LIST)):
			# get all hits
			raw_hits = [(n.start(0), n.end(0)) for n in re.finditer(KMER_LIST[ki], my_telseq)]
			coord_hit_dict.append({})
			for kmer_span in raw_hits:
				for j in range(kmer_span[0],kmer_span[1]):
					coord_hit_dict[-1][j] = True
					coord_hit_all[j]      = 1
			kmer_hit_dat[-1][0].append([n for n in raw_hits])
		#
		# remove hits of kmers that overlap with a hit from any of their superstrings
		#
		for ki in range(len(KMER_LIST)):
			del_list = []
			for ksi in range(len(kmer_hit_dat[-1][0][ki])):
				kmer_span = kmer_hit_dat[-1][0][ki][ksi]
				are_we_hit = False
				for sub_i in KMER_ISSUBSTRING[ki]:
					for j in range(kmer_span[0], kmer_span[1]):
						if j in coord_hit_dict[sub_i]:
							are_we_hit = True
							break
					if are_we_hit:
						break
				if are_we_hit:
					del_list.append(ksi)
			before_del_len = len(kmer_hit_dat[-1][0][ki])
			for di in sorted(del_list, reverse=True):
				del kmer_hit_dat[-1][0][ki][di]
			#print(ki, KMER_LIST[ki], KMER_ISSUBSTRING[ki], before_del_len, '-->', len(kmer_hit_dat[-1][0][ki]))
		#
		# collapse adjacent hits into larger blocks (so we have less polygons to plot)
		#
		for ki in range(len(KMER_LIST)):
			collapsed_kmer_spans = [[n[0],n[1]] for n in kmer_hit_dat[-1][0][ki]]
			for j in range(len(collapsed_kmer_spans)-1,0,-1):
				if collapsed_kmer_spans[j-1][1] == collapsed_kmer_spans[j][0]:
					collapsed_kmer_spans[j-1][1] = collapsed_kmer_spans[j][1]
					del collapsed_kmer_spans[j]
			kmer_hit_dat[-1][0][ki] = [n for n in collapsed_kmer_spans]
			#print(kmer_hit_dat[-1][0][ki])
	#
	plot_kmer_hits(kmer_hit_dat, KMER_COLORS, my_chr, OUT_PLOT)
	#

if __name__ == '__main__':
	main()
