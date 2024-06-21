# SeqRefFilter

When we try to filter some sequencing data by aligning it to a reference genome, sometimes we only want to use the ratio of matched bases to read length (match ratio) to determine whether to retain the read. This tool can help you calculate the distribution of match ratios in a BAM file, and after you decide on the threshold, it returns FASTQ file(s) from the BAM file.

Limitations:

1.	The degree of difference between sequencing data and the reference genome is uneven, and a single threshold may cause errors;
2.	Using the match ratio based on the entire read length may lead to the elimination of reads aligned to the edges of reference genome contigs or reads involving structural variations;
3. The tool is not suitable for filtering reads with low mapping quality.
4. The tool is suitable for short reads, such as Illumina reads.

## Installation

Recommend using bwa and samtools to align reads to the reference genome and convert the SAM file to the BAM file.

```bash
pip install seqreffilter
```

## Usage
mappping: Align reads to the reference genome.
```bash
bwa index reference.genome.fasta
bwa mem -o bwa_map.sam -t 8 reference.genome.fasta input_1_fq.gz input_2_fq.gz
samtools view -@ 8 -bS -o bwa_map.bam bwa_map.sam
```

stats: Calculate the distribution of match ratios in a BAM file.
```bash
SeqRefFilter stats -1 input_1_fq.gz -2 input_2_fq.gz -o stats bwa_map.bam
```

filter: Filter reads based on the match ratio.
```bash
SeqRefFilter filter -t 0.9 -o SeqRefFilter -1 input_1_fq.gz -2 input_2_fq.gz bwa_map.bam
```
