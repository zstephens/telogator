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

`python3 grab_reads_from_subtel_aln.py \ `  
`    --bam aln.bam \ `  
`    --fa clr-reads.fa.gz \ `  
`    --out subtel-reads.fa.gz \ `  
`    --bed resources/subtel_regions.bed `  

## (4 - PacBio CLR) align subtel reads to telogator reference:

`gunzip resources/t2t-telogator-ref.fa.gz`  

`pbmm2 align t2t-telogator-ref.fa subtel-reads.fa.gz subtel_aln.bam --preset SUBREAD --sort`  

## (4 - PacBio HiFi) align subtel reads to telogator reference:

We recommend using the [winnowmap](https://github.com/marbl/Winnowmap) aligner:

`winnowmap -W resources/repetitive_k15.txt -ax map-pb -Y t2t-telogator-ref.fa subtel-reads.fa.gz | samtools view -bh > subtel_aln-unsort.bam`  
`samtools sort -o subtel_aln.bam subtel_aln-unsort.bam`  

## (5) run telogator on subtel-only alignment:

`samtools view subtel_aln.bam | python3 telogator.py -i - -o telogator_out/`  

For Nanopore reads which may contain systematic errors in telomere regions, using the `-k` input option to select `resources/nanopore_kmers.tsv` may be advised. See [this preprint](https://www.biorxiv.org/content/10.1101/2022.01.11.475254v1) for more information.

## (6) create plots and output report:

`python3 merge_jobs.py -i telogator_out/`  

## Output files:

* `results.tsv`: telomere positions and lengths. column format:  
 * chr, boundary position, consensus TL, TL for each read, readlen for each read
* `tel_lens_violin_*.png`: violin plots of telomere lengths  
 * `all`: all reads
 * `chr`: only reads anchored to t2t reference
 * `alt`: only reads anchored to alternate assemblies

Additionally, the `--extra-readlen-plots` and `-rl` parameters can be used with `merge_jobs.py` to produce violin plots of read lengths, for QC purposes.
