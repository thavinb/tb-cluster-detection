#!/usr/bin/env python

'''
Convert tabular format of variants from
`gatk VariantsToTable` with `--genotype-fields GT`option
to FASTA format fo Multiple sequences alignment.
'''

import csv
import argparse
import sys


def vcftab2msa(table, output):
    with open(table, "r") as infile:
        tsv = csv.reader(infile, delimiter="\t")
        transposed_tsv = zip(*tsv)

    with open(output, "w") as outfile:
        # skip POS header
        transposed_tsv.__next__()
        for line in transposed_tsv:
            line = list(line)
            line[0] = ">" + line[0]
            line[0] = line[0].replace(".GT","\r\n")
            line.append("\r\n")
            text="".join(line)
            outfile.write(text)

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Turn gatk VariantsToTable tabular output format in to msa file",
        epilog="Example: python variantstabular_to_msa.py",
    )

    parser.add_argument(
        'filename',
        help="Genotypes GT tsv file from gatk VaraintsToTable"
    )

    parser.add_argument(
        'output',
        help="Multiple sequence alignment fasta"
    )

    return parser.parse_args(argv)

def main(argv=None):
    """
    Main function
    """
    args = parse_args(argv)
    vcftab2msa(args.filename, args.output)

if __name__ == "__main__":
    sys.exit(main())


