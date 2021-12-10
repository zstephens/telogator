import argparse
import os
import sys
import numpy as np

import matplotlib.pyplot as mpl

# parent directory of this script
telogator_base_dir = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
sys.path.append(telogator_base_dir)

from source.tg_util import exists_and_is_nonzero

def mad(x):
	med = np.median(x)
	return np.median(np.abs(x-med))

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23}

sorted_refs = [[n[1]+'_p', n[1]+'_q'] for n in sorted([(LEXICO_2_IND[k], k) for k in LEXICO_2_IND.keys()])]
sorted_refs = [val for sublist in sorted_refs for val in sublist]
sorted_refs = [n.replace('_','') for n in sorted_refs]

def main(raw_args=None):
	parser = argparse.ArgumentParser(description='simulate_reads.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
	parser.add_argument('-i', type=str, required=True,  metavar='input/',         help="* Input directory")
	parser.add_argument('-o', type=str, required=True,  metavar='accuracies.tsv', help="* Output simulation accuracies")
	args = parser.parse_args()

	OUT_BASE = args.i
	OUT_TSV  = args.o

	runs = sorted([n for n in os.listdir(OUT_BASE) if n[:4] == 'sim_'])

	ALL_COV = []
	ALL_LEN = []
	ALL_ERR = []
	for rd in runs:
		splt = rd.split('_')[1:]
		for n in splt:
			splt2 = n.split('-')
			if splt2[0] == 'c':
				ALL_COV.append(int(splt2[1]))
			elif splt2[0] == 'r':
				ALL_LEN.append(int(splt2[1]))
			elif splt2[0] == 'a':
				ALL_ERR.append(int(splt2[1]))
	ALL_COV = sorted(list(set(ALL_COV)))
	ALL_LEN = sorted(list(set(ALL_LEN)))
	ALL_ERR = sorted(list(set(ALL_ERR)))
	print('sweep:', ALL_COV, ALL_LEN, ALL_ERR)
	#
	ALL_TELPOS_DIST = {}
	ALL_TELLEN_DIST = {}
	ALL_FOUND_COUNT = {}
	for c in ALL_COV:
		ALL_TELPOS_DIST[c] = {}
		ALL_TELLEN_DIST[c] = {}
		ALL_FOUND_COUNT[c] = {}
		for r in ALL_LEN:
			ALL_TELPOS_DIST[c][r] = {}
			ALL_TELLEN_DIST[c][r] = {}
			ALL_FOUND_COUNT[c][r] = {}
			for a in ALL_ERR:
				ALL_TELPOS_DIST[c][r][a] = []
				ALL_TELLEN_DIST[c][r][a] = []
				ALL_FOUND_COUNT[c][r][a] = []

	for rd in runs:
		#
		truth_tlens = OUT_BASE + rd + '/tlen.tsv'
		if exists_and_is_nonzero(truth_tlens) == False:
			print(rd, 'missing tlen.tsv...')
			continue
		truth_dat = {}	# [chr] = (anchor_pos, tel_len)
		f = open(truth_tlens, 'r')
		for line in f:
			splt = line.strip().split('\t')
			truth_dat[splt[0].replace('_','')] = (int(splt[1]), int(splt[2]))
		f.close()
		#
		telogator_tlens = OUT_BASE + rd + '/telogator/results.tsv'
		if exists_and_is_nonzero(telogator_tlens) == False:
			print(rd, 'missing telogator/results.tsv...')
			continue
		telogator_dat = {}	# [chr] = (anchor_pos, [list, of, tel, lens, ...])
		f = open(telogator_tlens, 'r')
		for line in f:
			if line[0] != '#':
				splt = line.strip().split('\t')
				if splt[0][:3] == 'chr':
					tg_ref = splt[0].replace('_','')
				elif splt[0][:3] == 'alt':
					tg_ref = ''.join(splt[0].split('_')[:-1])
				tg_pos = int(splt[1])
				tg_len_list = [int(n) for n in splt[2].split(',')]
				if tg_ref not in telogator_dat:
					telogator_dat[tg_ref] = (tg_pos, tg_len_list)
				elif len(tg_len_list) > len(telogator_dat[tg_ref][1]):
					telogator_dat[tg_ref] = (tg_pos, tg_len_list)
				#if tg_ref not in telogator_dat:
				#	telogator_dat[tg_ref] = []
				#telogator_dat[tg_ref].append((int(splt[1]), [int(n) for n in splt[2].split(',')]))
		f.close()

		print('computing accuracy of run:', rd)
		pos_dists  = []
		len_dists  = []
		found_list = []
		for k in sorted_refs:
			gtp = truth_dat[k][0]
			gtl = truth_dat[k][1]
			if k in telogator_dat:
				#print(k, gtp, '~', telogator_dat[k][0], '=', abs(truth_dat[k][0]-telogator_dat[k][0]))
				mean_tlen = np.mean(telogator_dat[k][1])
				max_tlen  = np.max(telogator_dat[k][1])
				p90_tlen  = np.percentile(telogator_dat[k][1], 90)
				#print('---', gtl, ':', max_tlen, '=', int(abs(gtl - max_tlen)), ', reads =', len(telogator_dat[k][1]))
				#
				pos_dists.append(abs(truth_dat[k][0]-telogator_dat[k][0]))
				len_dists.append(int(abs(gtl - max_tlen)))
				found_list.append(1)
			elif 'alt'+k[3:] in telogator_dat:
				#print(k, 'FOUND IN ALTS')
				max_tlen = np.max(telogator_dat['alt'+k[3:]][1])
				len_dists.append(int(abs(gtl - max_tlen)))
				found_list.append(1)
			else:
				found_list.append(0)
		#
		my_key = [None, None, None]
		splt = rd.split('_')[1:]
		for n in splt:
			splt2 = n.split('-')
			if splt2[0] == 'c':
				my_key[0] = int(splt2[1])
			elif splt2[0] == 'r':
				my_key[1] = int(splt2[1])
			elif splt2[0] == 'a':
				my_key[2] = int(splt2[1])
		ALL_TELPOS_DIST[my_key[0]][my_key[1]][my_key[2]].extend([n for n in pos_dists])
		ALL_TELLEN_DIST[my_key[0]][my_key[1]][my_key[2]].extend([n for n in len_dists])
		ALL_FOUND_COUNT[my_key[0]][my_key[1]][my_key[2]].extend([n for n in found_list])
		#break
	print()

	f = open(OUT_TSV, 'w')
	header = 'cov' + '\t' +'readlen' + '\t' + 'error' + '\t' + 'mean_pos-error' + '\t' + 'mean_tlen-error' + '\t' + 'mean_found-frac'
	print(header)
	f.write(header + '\n')
	for c in ALL_COV:
		for r in ALL_LEN:
			for a in ALL_ERR:
				d1 = int(np.mean(ALL_TELPOS_DIST[c][r][a]))
				d2 = int(np.mean(ALL_TELLEN_DIST[c][r][a]))
				d3 = np.mean(ALL_FOUND_COUNT[c][r][a])
				print('\t'.join([str(n) for n in [c,r,a,d1,d2,d3]]))
				f.write('\t'.join([str(n) for n in [c,r,a,d1,d2,d3]]) + '\n')
	print()
	f.close()

if __name__ == '__main__':
	main()
