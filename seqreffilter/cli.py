import argparse
from seqreffilter.src.pipeline import stats_main, filter_main


class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='SeqRefFilter',
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for stats
        parser_a = subparsers.add_parser('stats',
                                         description='calculate the distribution of match ratios in a BAM file\n')

        parser_a.add_argument('input_bam', type=str,
                                help='input bam file')
        parser_a.add_argument('-1', '--input_fq_1', type=str,
                                help='input fastq file 1')
        parser_a.add_argument('-2', '--input_fq_2', type=str,
                                help='input fastq file 2')
        parser_a.add_argument('-o', '--output', type=str,
                                help='output prefix', default='stats')
        parser_a.set_defaults(func=stats_main)

        # argparse for filter
        parser_b = subparsers.add_parser('filter',
                                         description='filter the reads in a BAM file\n')
        
        parser_b.add_argument('input_bam', type=str,
                                help='input bam file')
        parser_b.add_argument('-1', '--input_fq_1', type=str,
                                help='input fastq file 1')
        parser_b.add_argument('-2', '--input_fq_2', type=str,
                                help='input fastq file 2')
        parser_b.add_argument('-o', '--output', type=str,
                                help='output prefix', default='filtered')
        parser_b.add_argument('-t', '--threshold', type=float, default=0.8,
                                help='threshold of match ratio, default is 0.8')
        parser_b.set_defaults(func=filter_main)

        self.arg_parser = parser

        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()

        if self.args.subcommand_name == 'stats':
            stats_main(self.args)
        elif self.args.subcommand_name == 'filter':
            filter_main(self.args)
        else:
            self.arg_parser.print_help()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
