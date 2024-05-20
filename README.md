# Telogator --> Telogator2

This software is now deprecated. Please use [Telogator2](https://github.com/zstephens/telogator2).

# Telogator
A method for measuring chromosome-specific telomere length (TL) from long reads.

If this software has been useful for your work, please cite us at:

Stephens, Zachary, et al. "Telogator: a method for reporting chromosome-specific telomere lengths from long reads." *Bioinformatics* 38.7 (2022): 1788-1793. https://doi.org/10.1093/bioinformatics/btac005


## Dependencies:

Telogator dependencies can be easily installed via conda:

```bash
# create conda environment
conda env create -f conda_env_telogator.yaml

# activate environment
conda activate telogator
```

## (1) making T2T + alt subtel reference:

Obtain t2t reference sequence from https://github.com/marbl/CHM13 (`chm13v2.0.fa`)

Append alternate subtelomere assemblies:

```bash
cat chm13v2.0.fa resources/stong_subtels.fa > t2t-and-subtel.fa
samtools faidx t2t-and-subtel.fa
```


## (2) align reads to whole genome reference:

For PacBio CLR reads:

```bash
pbmm2 align t2t-and-subtel.fa clr-reads.fa.gz aln.bam --preset SUBREAD --sort
```

For PacBio HiFi reads:

```bash
pbmm2 align t2t-and-subtel.fa hifi-reads.fa.gz aln.bam --preset HiFi --sort
```

For Oxford Nanopore reads:

```bash
minimap2 -ax map-ont -t 6 -N 5 -Y -L -o aln-unsort.sam t2t-and-subtel.fa ont-reads.fa.gz
samtools sort -o aln.bam aln-unsort.sam
samtools index aln.bam
```


## (3) extract subset of reads that were mapped to subtelomere regions:

```bash
python3 get_subtel_reads.py -b aln.bam -i clr-reads.fa.gz -o subtel-reads.fa.gz
```


## (4) align read subset to telogator reference:

We recommend using the [winnowmap](https://github.com/marbl/Winnowmap) aligner for best results. An example for HiFi reads:

```bash
winnowmap -ax map-hifi -Y -W resources/repetitive-k15.txt -o subtel_aln-unsort.sam t2t-telogator-ref.fa subtel-reads.fa.gz
samtools sort -o subtel_aln.bam subtel_aln-unsort.sam
samtools index subtel_aln.bam
```

Though this alignment could be done using the same tools as in step 2, if desired.



## (5) run telogator on subtel-only alignment:

```bash
python3 telogator.py -i subtel_aln.bam -o telogator_out/
```

The `-p` input option specifies the number of processes to use (default: 4).

For Nanopore reads, which might have systematic sequencing errors in telomere regions, try adding the `-r ont` input option if Telogator is not finding many telomeres.

