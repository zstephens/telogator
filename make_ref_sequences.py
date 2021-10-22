import os
import sys
import time
import numpy as np

from source.tg_util import makedir, rm, exists_and_is_nonzero, read_fq_entry, index_ref, read_ref_py3, RC

LEXICO_2_IND = {'chr1':1, 'chr2':2, 'chr3':3, 'chr10':10, 'chr11':11, 'chr12':12, 'chr19':19, 'chr20':20,
				'chr4':4, 'chr5':5, 'chr6':6, 'chr13':13, 'chr14':14, 'chr15':15, 'chr21':21, 'chr22':22,
				'chr7':7, 'chr8':8, 'chr9':9, 'chr16':16, 'chr17':17, 'chr18':18, 'chrX' :23, 'chrY' :24, 'chrM' :25}

STONG_CONTIGS = [('1ptel_1-500K_1_12_12v2',     'alt1_p_0'),
                 ('1qtel_1-500K_1_12_12',       'alt1_q_0'),
                 ('2ptel_1-500K_1_12_12',       'alt2_p_0'),
                 ('2qtel_1-500K_1_12_12',       'alt2_q_0'),
                 ('3ptel_1-500K_1_12_12',       'alt3_p_0'),
                 ('3qtel_1-500K_1_12_12',       'alt3_q_0'),
                 ('4ptel_1-500K_1_12_12',       'alt4_p_0'),
                 ('4qtel_1-500K_1_12_12',       'alt4_q_0'),
                 ('5ptel_1-500K_1_12_12',       'alt5_p_0'),
                 ('5qtel_1-500K_1_12_12',       'alt5_q_0'),
                 ('6ptel_1-500K_1_12_12',       'alt6_p_0'),
                 ('6qtel_1-500K_1_12_12',       'alt6_q_0'),
                 ('7ptel_short_1-500K_1_12_12', 'alt7_p_0'),
                 ('7qtel_1-500K_1_12_12',       'alt7_q_0'),
                 ('8ptel_1_500K_1_12_12',       'alt8_p_0'),
                 ('8qtel_short_1-500K_1_12_12', 'alt8_q_0'),
                 ('9ptel_1-500K_1_12_12',       'alt9_p_0'),
                 ('9qtel_1-500K_1_12_12',       'alt9_q_0'),
                 ('10ptel_1-500K_1_12_12',      'alt10_p_0'),
                 ('10qtel_1-500K_1_12_12',      'alt10_q_0'),
                 ('11ptel_1-500K_1_12_12',      'alt11_p_0'),
                 ('11qtel_1-500K_1_12_12',      'alt11_q_0'),
                 ('12ptel_1-500K_1_12_12',      'alt12_p_0'),
                 ('12qtel_1-500K_1_12_12',      'alt12_q_0'),
                 ('13qtel_1-500K_new_4_3_12_hg19orientation', 'alt13_q_0'),
                 ('14qtel_1-500K_1_12_12',                    'alt14_q_0'),
                 ('15qtel_1-500K_1_12_12',                    'alt15_q_0'),
                 ('16ptel_short_1-500K_1_12_12',      'alt16_p_0'),
                 ('16qtel_1-500K_1_12_12',            'alt16_q_0'),
                 ('17ptel_1_500K_1_12_12',            'alt17_p_0'),
                 ('17qtel_1-500K_1_12_12v2',          'alt17_q_0'),
                 ('18ptel_1-500K_1_12_12',            'alt18_p_0'),
                 ('18qtel_1-500K_1_12_12',            'alt18_q_0'),
                 ('19ptel_1-500K_1_12_12',            'alt19_p_0'),
                 ('19qtel_1-500K_1_12_12',            'alt19_q_0'),
                 ('20ptel_1-500K_1_12_12',            'alt20_p_0'),
                 ('20qtel_1-500Kshortallele_1_12_12', 'alt20_q_0'),
                 ('21qtel_1-500K_1_12_12',  'alt21_q_0'),
                 ('22qtel_1-500_1_12_12',   'alt22_q_0'),
                 ('Xqtel_1-500K_1_12_12',   'altX_q_0'),
                 ('Yqtel_1-500K_1_12_12',   'altY_q_0'),
                 ('XpYptel_1-500K_1_12_12', 'altX_p_0')]



rfn = '/Users/zach/Downloads/Supplemental_FileS1.txt'
ref_inds = index_ref(rfn)
ri_rev   = {ref_inds[n][0]:n for n in range(len(ref_inds))}

FASTA_WIDTH = 60
SUBTEL_BUFF = 500000
ofn = '/Users/zach/Desktop/stong_subtel_new.fa'
f   = open(ofn, 'w')
for i in range(len(STONG_CONTIGS)):
	print(STONG_CONTIGS[i][0], ri_rev[STONG_CONTIGS[i][0]])
	my_ri     = ri_rev[STONG_CONTIGS[i][0]]
	ref_seq   = read_ref_py3(rfn, ref_inds[my_ri])
	f.write('>' + STONG_CONTIGS[i][1] + '\n')
	for j in range(0,len(ref_seq),FASTA_WIDTH):
		f.write(ref_seq[j:j+FASTA_WIDTH] + '\n')
f.close()

rfn = '/Users/zach/Desktop/mayo/ref/chm13.draft_v1.1.fa'
bfn = '/Users/zach/Downloads/chm13.draft_v1.1.telomere.bed'
ofn = '/Users/zach/Desktop/t2t_subtels_maskedtel.fa'

out_bed_new = '/Users/zach/Desktop/subtel_regions_t2t_plus_alts.bed'

ref_inds = index_ref(rfn)
ri_rev   = {ref_inds[n][0]:n for n in range(len(ref_inds))}

bed_regions = {}
f = open(bfn, 'r')
for line in f:
	splt = line.strip().split('\t')
	if splt[0] not in bed_regions:
		bed_regions[splt[0]] = []
	bed_regions[splt[0]].append((int(splt[1]), int(splt[2])))
f.close()
sorted_refs = [n[1] for n in sorted([(LEXICO_2_IND[k],k) for k in bed_regions.keys()])]

f  = open(ofn, 'w')
f2 = open(out_bed_new, 'w')
for ref in sorted_refs:
	my_ri     = ri_rev[ref]
	ref_seq   = read_ref_py3(rfn, ref_inds[my_ri])
	#
	for i in range(len(bed_regions[ref])):
		n = bed_regions[ref][i]
		n0 = max([0, n[0]])
		n1 = min([len(ref_seq), n[1]])
		ref_seq = ref_seq[:n0] + 'N'*(n1-n0) + ref_seq[n1:]
	#
	ref_seq_p = ref_seq[:SUBTEL_BUFF]
	ref_seq_q = ref_seq[len(ref_seq)-SUBTEL_BUFF:]
	f2.write(ref + '\t' + '1' + '\t' + str(SUBTEL_BUFF) + '\t' + ref + '_p' + '\n')
	f2.write(ref + '\t' + str(len(ref_seq)-SUBTEL_BUFF) + '\t' + str(len(ref_seq)) + '\t' + ref + '_q' + '\n')
	#
	my_ref_name = ref + '_p'
	f.write('>' + my_ref_name + '\n')
	for i in range(0,len(ref_seq_p),FASTA_WIDTH):
		f.write(ref_seq_p[i:i+FASTA_WIDTH] + '\n')
	#
	my_ref_name = ref + '_q'
	f.write('>' + my_ref_name + '\n')
	for i in range(0,len(ref_seq_q),FASTA_WIDTH):
		f.write(ref_seq_q[i:i+FASTA_WIDTH] + '\n')
#
for i in range(len(STONG_CONTIGS)):
	f2.write(STONG_CONTIGS[i][1] + '\t' + '1' + '\t' + '500000' + '\n')
f2.close()
f.close()

