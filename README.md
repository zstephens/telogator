# Telogator
A method for measuring chromosome-specific TL from long reads.

If this software has been useful for your work, please cite us at:

Stephens, Zachary, et al. "Telogator: a method for reporting chromosome-specific telomere lengths from long reads." *Bioinformatics* 38.7 (2022): 1788-1793. https://doi.org/10.1093/bioinformatics/btac005

## (1) making T2T + alt subtel reference:

Obtain t2t reference sequence from https://github.com/marbl/CHM13 (`chm13v2.0.fa`)

Append alternate subtelomere assemblies:

`cat chm13v2.0.fa resources/stong_subtels.fa > t2t-and-subtel.fa`  
`samtools faidx t2t-and-subtel.fa`  

## (2) align reads to whole genome reference:

For PacBio CLR reads:

`pbmm2 align t2t-and-subtel.fa clr-reads.fa.gz aln.bam --preset SUBREAD --sort`  

For PacBio HiFi reads:

`pbmm2 align t2t-and-subtel.fa hifi-reads.fa.gz aln.bam --preset HiFi --sort`  

For Oxford Nanopore long (noisy) reads:

`minimap2 -ax map-ont -t 6 -N 5 -Y -L t2t-and-subtel.fa ont-reads.fa.gz | samtools view -bhS - > aln-unsort.bam`  
`samtools sort -o aln.bam aln-unsort.bam`  
`samtools index aln.bam`  

## (3) get subset of reads that were mapped to subtelomere regions:

`python3 get_subtel_reads.py \ `  
`    --bam aln.bam \ `  
`    --in-reads clr-reads.fa.gz \ `  
`    --out-reads subtel-reads.fa.gz`  

## (4) align read subset to telogator reference:

We recommend using the [winnowmap](https://github.com/marbl/Winnowmap) aligner for best results:

`winnowmap -W resources/repetitive_k15.txt \ `  
`    -ax map-pb \ `  
`    -Y \ `  
`    t2t-telogator-ref.fa \ `  
`    subtel-reads.fa.gz \ `  
`    | samtools view -bh > subtel_aln-unsort.bam`  
`samtools sort -o subtel_aln.bam subtel_aln-unsort.bam`  
`samtools index subtel_aln.bam`  

Though this step could be done using the same aligners as in step 2, if desired.

## (5) run telogator on subtel-only alignment:

`python3 telogator.py -i subtel_aln.bam -o telogator_out/`  

The `-p` input option specifies the number of processes to use (default: 4).
