import pywt
import re
import regex

import numpy as np

from source.tg_prob import DiscreteDistribution
from source.tg_util import RC

COEF_EDIT_0 = 2.
COEF_EDIT_1 = 1.

#
#
def get_telomere_kmer_density(read_dat, kmer_list, tel_window, smoothing=False):
	re_hits = [[], []]
	for i in range(len(kmer_list)):
		re_hits[0].extend([(n.start(0), n.end(0), i, 0) for n in re.finditer(kmer_list[i], read_dat)])
		re_hits[1].extend([(n.start(0), n.end(0), i, 1) for n in regex.finditer("("+kmer_list[i]+"){e<=1}", read_dat, overlapped=True)])
	#
	tel_hit_p0 = np.zeros(len(read_dat))
	tel_hit_p1 = np.zeros(len(read_dat))
	for re_hit in re_hits[0]:
		tel_hit_p0[re_hit[0]:re_hit[1]] = 1
	for re_hit in re_hits[1]:
		tel_hit_p1[re_hit[0]:re_hit[1]] = 1
	tel_hit_cum_e0 = np.cumsum(tel_hit_p0)
	tel_hit_cum_e1 = np.cumsum(tel_hit_p1)
	#
	tel_density_e0 = []
	tel_density_e1 = []
	for i in range(len(read_dat) - tel_window):
		tel_density_e0.append(float(tel_hit_cum_e0[i+tel_window] - tel_hit_cum_e0[i]) / tel_window)
		tel_density_e1.append(float(tel_hit_cum_e1[i+tel_window] - tel_hit_cum_e1[i]) / tel_window)
	#
	if smoothing:
		tel_density_e0 = wavelet_smooth(tel_density_e0)
		tel_density_e1 = wavelet_smooth(tel_density_e1)
	#
	return (tel_density_e0, tel_density_e1)

#
#
def get_telomere_regions(td_p_e0, td_p_e1, td_q_e0, td_q_e1, tel_window, pthresh, smoothing=True):
	min_win_pos  = min([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
	max_win_pos  = max([len(td_p_e0), len(td_p_e1), len(td_q_e0), len(td_q_e1)])
	p_vs_q_power = np.zeros(max_win_pos)
	tel_regions  = [[0,None,None]]
	#
	for i in range(min_win_pos):
		c0 = td_p_e0[i] - td_q_e0[i]
		c1 = td_p_e1[i] - td_q_e1[i]
		p_vs_q_power[i] = (COEF_EDIT_0*(c0) + COEF_EDIT_1*(c1)) / float(COEF_EDIT_0 + COEF_EDIT_1)
	if smoothing:
		p_vs_q_power = wavelet_smooth(p_vs_q_power)
	p_vs_q_power = p_vs_q_power[:-tel_window]
	#
	if p_vs_q_power[0] >= pthresh:
		tel_regions[-1][2] = 'p'
	elif p_vs_q_power[0] <= -pthresh:
		tel_regions[-1][2] = 'q'
	for i in range(1,len(p_vs_q_power)):
		if tel_regions[-1][2] == 'p' and p_vs_q_power[i] >= pthresh:
			tel_regions[-1][1] = i+1
		elif tel_regions[-1][2] == 'q' and p_vs_q_power[i] <= -pthresh:
			tel_regions[-1][1] = i+1
		elif tel_regions[-1][2] == None and abs(p_vs_q_power[i]) < pthresh:
			tel_regions[-1][1] = i+1
		else:
			if tel_regions[-1][1] == None:	# previous region was only 1bp wide, ugh
				tel_regions[-1][1] = tel_regions[-1][0]+1
			tel_regions.append([i,None,None])
			if p_vs_q_power[i] >= pthresh:
				tel_regions[-1][2] = 'p'
			elif p_vs_q_power[i] <= -pthresh:
				tel_regions[-1][2] = 'q'
	if tel_regions[-1][1] == None:
		tel_regions[-1][1] = len(p_vs_q_power)-1
	#
	return (p_vs_q_power, tel_regions)

#
# median absolute deviation, with some normalization
#
def mad(x, c=1.4826):
	return c * np.median(np.abs(x - np.median(x)))

#
# adapted from: http://connor-johnson.com/2016/01/24/using-pywavelets-to-remove-high-frequency-noise/
#
def wavelet_smooth(x, wavelet="db4", smoothing_level=5, fail_sigma=1e-3):
	# calculate the wavelet coefficients
	coeff = pywt.wavedec( x, wavelet, mode="per" )
	# calculate a threshold
	sigma = mad( coeff[-smoothing_level] )
	# if sigma below some value (i.e. signal is mostly 0s) return unsmoothed signal
	if sigma < fail_sigma:
		return x
	# changing this threshold also changes the behavior, but I have not played with this very much
	uthresh = sigma * np.sqrt( 2*np.log( len( x ) ) )
	coeff[1:] = ( pywt.threshold( i, value=uthresh, mode="soft" ) for i in coeff[1:] )
	# reconstruct the signal using the thresholded coefficients
	y = pywt.waverec( coeff, wavelet, mode="per" )
	return y

#
#
def write_kmer_prob_matrix(fn, kmer_list, trans_matrix):
	f = open(fn, 'w')
	for k in kmer_list:
		f.write('#\t' + k + '\n')
	for i in range(len(kmer_list)):
		f.write('\t'.join([str(n) for n in trans_matrix[i,:].tolist()]) + '\n')
	f.close()

#
#
def read_kmer_prob_matrix(fn, return_trans_dict=False):
	keys = []
	prob = []
	f = open(fn, 'r')
	for line in f:
		splt = line.strip().split('\t')
		if splt[0][0] == '#':
			keys.append(splt[1])
		else:
			prob.append([float(n) for n in splt])
	f.close()
	prob = np.array(prob)
	if return_trans_dict:
		trans_dict = {}
		for i in range(len(keys)):
			trans_dict[keys[i]] = DiscreteDistribution(prob[i,:].tolist(), keys)
		return trans_dict
	else:
		return (keys, prob)

#
#
def sample_telomere(trans_dist_dict, length, tel_type='p', init_string='TAACCC'):
	out_s = init_string
	current_kmer = init_string
	while len(out_s) < length:
		next_kmer = trans_dist_dict[current_kmer].sample()
		out_s += next_kmer
		current_kmer = next_kmer
	if tel_type == 'p':
		return out_s
	elif tel_type == 'q':
		return RC(out_s)
	else:
		print('Error: unknown tel type in sample_telomere()')
		exit(1)

#
#
def get_nonoverlapping_kmer_hits(my_telseq, KMER_LIST, KMER_ISSUBSTRING):
	out_dat = []
	coord_hit_dict = []
	for kmer_list_i in range(len(KMER_LIST)):
		# get all hits
		raw_hits = [(n.start(0), n.end(0)) for n in re.finditer(KMER_LIST[kmer_list_i], my_telseq)]
		coord_hit_dict.append({})
		for kmer_span in raw_hits:
			for j in range(kmer_span[0],kmer_span[1]):
				coord_hit_dict[-1][j] = True
		out_dat.append([n for n in raw_hits])
	#
	# remove hits of kmers that overlap with a hit from any of their superstrings
	#
	for kmer_list_i in range(len(KMER_LIST)):
		del_list = []
		for ksi in range(len(out_dat[kmer_list_i])):
			kmer_span = out_dat[kmer_list_i][ksi]
			are_we_hit = False
			for sub_i in KMER_ISSUBSTRING[kmer_list_i]:
				for j in range(kmer_span[0], kmer_span[1]):
					if j in coord_hit_dict[sub_i]:
						are_we_hit = True
						break
				if are_we_hit:
					break
			if are_we_hit:
				del_list.append(ksi)
		before_del_len = len(out_dat[kmer_list_i])
		for di in sorted(del_list, reverse=True):
			del out_dat[kmer_list_i][di]
	#
	# collapse adjacent hits into larger blocks (so we have less polygons to plot)
	#
	for kmer_list_i in range(len(KMER_LIST)):
		collapsed_kmer_spans = [[n[0],n[1]] for n in out_dat[kmer_list_i]]
		for j in range(len(collapsed_kmer_spans)-1,0,-1):
			if collapsed_kmer_spans[j-1][1] == collapsed_kmer_spans[j][0]:
				collapsed_kmer_spans[j-1][1] = collapsed_kmer_spans[j][1]
				del collapsed_kmer_spans[j]
		out_dat[kmer_list_i] = [n for n in collapsed_kmer_spans]
		#print(out_dat[kmer_list_i])
	#
	return out_dat
