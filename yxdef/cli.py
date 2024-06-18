import argparse
from yxdef.src.pipeline import subcmd1_main, subcmd2_main


class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='yxdef',
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for subcmd1
        parser_a = subparsers.add_parser('subcmd1',
                                         description='Description for subcmd1\n')

        parser_a.add_argument('input', type=str,
                                help='input file')
        parser_a.add_argument('output', type=str,
                                help='output file')
        parser_a.add_argument('-t', '--threads', type=int,
                                help='number of threads, default 1', default=1)
        parser_a.add_argument('-d', '--dry_run', action='store_true',
                                help='dry run, default False')
        parser_a.set_defaults(func=subcmd1_main)

        # argparse for subcmd2
        parser_b = subparsers.add_parser('subcmd2',
                                         description='Description for subcmd2\n')
        
        parser_b.add_argument('input', type=str,
                                help='input file')
        parser_b.add_argument('output', type=str,
                                help='output file')
        parser_b.add_argument('-t', '--threads', type=int,
                                help='number of threads, default 1', default=1)
        parser_b.add_argument('-d', '--dry_run', action='store_true',
                                help='dry run, default False')
        parser_b.set_defaults(func=subcmd2_main)

        self.arg_parser = parser

        self.args = parser.parse_args()

    def run(self):
        self.run_arg_parser()

        if self.args.subcommand_name == 'subcmd1':
            subcmd1_main(self.args)
        elif self.args.subcommand_name == 'subcmd2':
            subcmd2_main(self.args)
        else:
            self.arg_parser.print_help()


def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()
