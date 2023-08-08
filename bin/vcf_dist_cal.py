#!/usr/bin/env python

'''
vcf_dist_cal.py <vcf1> <vcf2>
Calculate SNP distant base on ALT of vcfs
This expect vcf output from gatk HaploTypeCaller
run with following setting
--output-mode EMIT_VARIANTS_ONLY
-ploidy 1
'''

import vcf
from vcf import utils
import gzip
import sys
import argparse

def cal_dist(vcf1, vcf2):
    """
    vcf.utils.walk_together give the list of vcf parsed data on each position
    [ None      , <vcf2.obj>] SNP on vcf1 is REF and vcf2 is ALT at this POS
    [ <vcf1.obj>, None      ] SNP on vcf1 is ALT and vcf2 is REF at this POS
    [ <vcf1.obj>, <vcf2.obj>] SNP on vcf1 is ALT and vcf2 is ALT at this POS

    iterate over each record (POS) and compare base
    """
    vcf1_Reader = vcf.Reader(filename=vcf1)
    vcf2_Reader = vcf.Reader(filename=vcf2)
    conn_vcf    = vcf.utils.walk_together(vcf1_Reader, vcf2_Reader)


    dist = 0
    for row in conn_vcf:
        # None mean same base as REF
        pos = row[0].POS if row[0] != None else row[1].POS
        v1  = row[0].ALT if row[0] != None else row[1].REF
        v2  = row[1].ALT if row[1] != None else row[0].REF
        if v1 != v2:
            dist += 1
            # print(pos, v1, v2)
    print(dist)



def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Calculate SNP distant base on ALT of vcfs",
        epilog="Example: vcf_dist_cal.py <vcf1> <vcf2>",
    )

    parser.add_argument(
        'vcf1',
        help="VCF file 1"
    )

    parser.add_argument(
        'vcf2',
        help="VCF file 1"
    )

    return parser.parse_args(argv)

def main(argv=None):
    """
    Main function
    """
    args = parse_args(argv)
    cal_dist(args.vcf1, args.vcf2)

if __name__ == "__main__":
    sys.exit(main())


