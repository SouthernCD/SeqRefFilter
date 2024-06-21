from seqreffilter.src.utils import parse_bam_by_reads
import numpy as np
from yxseq import read_fastq_big
from yxmath.plot import hist


def stats_main(args):

    if args.input_fq_1.endswith('.gz'):
        fastq_reader = read_fastq_big(args.input_fq_1, args.input_fq_2, True)
    else:
        fastq_reader = read_fastq_big(args.input_fq_1, args.input_fq_2, False)

    # get the match ratios
    match_ratios = []
    for map_list in parse_bam_by_reads(args.input_bam):
        r_name = map_list[0].query_name
        r = next(fastq_reader)

        if r_name != r[0].seqname:
            raise ValueError(
                "Read name in bam file and fastq file are not match")

        r1_len = len(r[0].seq)
        if len(r) > 1:
            r2_len = len(r[1].seq)
        else:
            r2_len = 0

        best_r1_match = 0
        best_r2_match = 0

        for m in map_list:

            if m.is_read1:
                if m.is_unmapped:
                    f_match = 0
                else:
                    f_match = m.get_cigar_stats()[0][0]
                if f_match > best_r1_match:
                    best_r1_match = f_match

            if m.is_read2:
                if m.is_unmapped:
                    f_match = 0
                else:
                    f_match = m.get_cigar_stats()[0][0]
                if f_match > best_r2_match:
                    best_r2_match = f_match

        match_ratio = (best_r1_match + best_r2_match) / (r1_len + r2_len)

        match_ratios.append(match_ratio)

    # statistics
    match_ratios = np.array(match_ratios)
    n, bins = hist(match_ratios, save_file=args.output + '.png')

    with open(args.output + '.txt', 'w') as f:
        for i in range(len(n)):
            f.write(f'{bins[i]}\t{bins[i+1]}\t{n[i]}\n')


def filter_main(args):
    if args.input_fq_1.endswith('.gz'):
        fastq_reader = read_fastq_big(args.input_fq_1, args.input_fq_2, True)
    else:
        fastq_reader = read_fastq_big(args.input_fq_1, args.input_fq_2, False)

    fq1 = open(args.output + '_1.fq', 'w')
    if args.input_fq_2:
        fq2 = open(args.output + '_2.fq', 'w')

    for map_list in parse_bam_by_reads(args.input_bam):
        r_name = map_list[0].query_name
        r = next(fastq_reader)

        if r_name != r[0].seqname:
            raise ValueError(
                "Read name in bam file and fastq file are not match")

        r1_len = len(r[0].seq)
        if len(r) > 1:
            r2_len = len(r[1].seq)
        else:
            r2_len = 0

        best_r1_match = 0
        best_r2_match = 0

        for m in map_list:

            if m.is_read1:
                if m.is_unmapped:
                    f_match = 0
                else:
                    f_match = m.get_cigar_stats()[0][0]
                if f_match > best_r1_match:
                    best_r1_match = f_match

            if m.is_read2:
                if m.is_unmapped:
                    f_match = 0
                else:
                    f_match = m.get_cigar_stats()[0][0]
                if f_match > best_r2_match:
                    best_r2_match = f_match

        match_ratio = (best_r1_match + best_r2_match) / (r1_len + r2_len)

        if match_ratio >= args.threshold:
            fq1.write(f'@{r[0].seqname}\n{r[0].seq}\n+\n{r[0].quality}\n')
            if len(r) > 1:
                fq2.write(f'@{r[1].seqname}\n{r[1].seq}\n+\n{r[1].quality}\n')

    fq1.close()
    if args.input_fq_2:
        fq2.close()


if __name__ == '__main__':
    pass