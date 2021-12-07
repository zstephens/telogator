# Telogator
A method for measuring chromosome-specific TL from long reads.

## (1) making T2T + alt subtel reference:

Obtain t2t reference sequence from https://github.com/marbl/CHM13 (`chm13.draft_v1.1.fa`)

Append alternate subtelomere assemblies:

`gunzip resources/stong_subtels.fa.gz`  

`cat chm13.draft_v1.1.fa stong_subtels.fa > t2t-and-subtel.fa`  

`samtools faidx t2t-and-subtel.fa`  

## (2) align reads to whole genome (PacBio CLR):

`pbmm2 align t2t-and-subtel.fa clr-reads.fa.gz aln.bam --preset SUBREAD --sort`  

## (3) extract reads mapped to subtels:

`python3 grab_subreads_from_t2t-and-subtel_aln.py \ `  
`    --bam aln.bam \ `  
`    --fa clr-reads.fa.gz \ `  
`    --out subtel-reads.fa.gz \ `  
`    --bed resources/subtel_regions.bed `  

## (4) align subtel reads to subtel-only reference (PacBio CLR):

`gunzip resources/t2t-telogator-ref.fa.gz`  

`pbmm2 align t2t-telogator-ref.fa subtel-reads.fa.gz subtel_aln.bam --preset SUBREAD --sort`  

## (5) run telogator on subtel-only alignment:

`samtools view subtel_aln.bam | python3 telogator.py -i - -o telogator_out/`  

## (6) create plots and output report:

`python3 merge_jobs.py -i telogator_out/`  
