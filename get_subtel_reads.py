#!/usr/bin/env python
import argparse
import gzip
import pathlib
import pysam

from source.tg_reader import TG_Reader
from source.tg_util import exists_and_is_nonzero


def get_clr_readname(rn):
    return '/'.join(rn.split('/')[:-1])


def get_sra_readname(rn):
    return rn.split(' ')[0]


def get_ccs_readname(rn):
    return rn.split(' ')[0]


def main(raw_args=None):
    parser = argparse.ArgumentParser(description='grab_subreads_from_t2t-and-subtel_aln.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    parser.add_argument('--bam',       type=str, required=True,  metavar='input.bam',        help="* Input BAM")
    parser.add_argument('--in-reads',  type=str, required=True,  metavar='input.fa',         help="* Input reads (.fa / .fa.gz / .fq / .fq.gz)")
    parser.add_argument('--out-reads', type=str, required=True,  metavar='output.fa',        help="* Output reads (.fa / .fa.gz / .fq / .fq.gz)")
    parser.add_argument('--bed',       type=str, required=False, metavar='subtel.bed',       help="Subtel regions",   default='')
    parser.add_argument('--readtype',  type=str, required=False, metavar='CCS / CLR / SRA',  help="Read name format", default='SRA')
    args = parser.parse_args()

    IN_BAM     = args.bam
    IN_READS   = args.in_reads
    OUT_READS  = args.out_reads
    #
    SUBTEL_BED = args.bed
    READTYPE   = args.readtype

    if SUBTEL_BED == '':
        print('using default subtelomere regions...')
        sim_path  = pathlib.Path(__file__).resolve().parent
        SUBTEL_BED = str(sim_path) + '/resources/subtel_regions.bed'

    if exists_and_is_nonzero(IN_BAM) is False:
        print('Error: input.bam not found.')
        exit(1)
    if exists_and_is_nonzero(IN_READS) is False:
        print('Error: input.fa not found.')
        exit(1)
    if exists_and_is_nonzero(SUBTEL_BED) is False:
        print('Error: subtel.bed not found.')
        exit(1)

    samfile = pysam.AlignmentFile(IN_BAM, "rb")

    OUTPUT_IS_FASTQ = False
    if OUT_READS[-3:].lower() == '.fq' or OUT_READS[-6:].lower() == '.fq.gz':
        OUTPUT_IS_FASTQ = True
    INPUT_IS_FASTQ = False
    if IN_READS[-3:].lower() == '.fq' or IN_READS[-6:].lower() == '.fq.gz':
        INPUT_IS_FASTQ = True
    if INPUT_IS_FASTQ is False and OUTPUT_IS_FASTQ is True:
        print('Error: input is fasta and output is fastq.')
        exit(1)

    bed_str = []
    f = open(SUBTEL_BED, 'r')
    for line in f:
        splt = line.strip().split('\t')
        bed_str.append([n for n in splt[:3]])
    f.close()

    OUT_DIR = '/'.join(OUT_READS.split('/')[:-1])
    if len(OUT_DIR) == 0:
        OUT_DIR = '.'
    OUT_DIR += '/'

    print('getting readnames from bam...')
    rn_dict = {}
    for i in range(len(bed_str)):
        bed_dat = bed_str[i][0] + ':' + bed_str[i][1] + '-' + bed_str[i][2]
        print('-', bed_dat)
        try:
            for aln in samfile.fetch(region=bed_dat):
                sam_line = str(aln).split('\t')
                my_rnm   = sam_line[0]
                if READTYPE == 'CLR':
                    rn_dict[get_clr_readname(my_rnm)] = True
                elif READTYPE == 'SRA':
                    rn_dict[get_sra_readname(my_rnm)] = True
                elif READTYPE == 'CCS':
                    rn_dict[get_ccs_readname(my_rnm)] = True
                else:
                    print('Error: unknown read type, must be: CCS / CLR / SRA')
                    exit(1)
        except ValueError:
            print('skipping contig that was not present in input BAM:', bed_str[i][0])

    if OUT_READS[-3:].lower() == '.gz':
        f_out = gzip.open(OUT_READS, 'wt')
    else:
        f_out = open(OUT_READS, 'w')
    #
    my_reader = TG_Reader(IN_READS)
    #
    while True:
        (my_name, my_rdat, my_qdat) = my_reader.get_next_read()
        if not my_name:
            break
        got_hits = [(READTYPE == 'CLR' and get_clr_readname(my_name) in rn_dict),
                    (READTYPE == 'SRA' and get_sra_readname(my_name) in rn_dict),
                    (READTYPE == 'CCS' and get_ccs_readname(my_name) in rn_dict)]
        if any(got_hits):
            if OUTPUT_IS_FASTQ:
                f_out.write('@' + my_name + '\n' + my_rdat + '\n' + '+' + '\n' + my_qdat + '\n')
            else:
                f_out.write('>' + my_name + '\n' + my_rdat + '\n')
    #
    my_reader.close()
    f_out.close()


if __name__ == '__main__':
    main()
