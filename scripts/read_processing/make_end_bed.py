#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import sys

QNAME=0
FLAG=1
RNAME=2
POS=3
TLEN=8
SEQ=9

#cutwin=10

parser = argparse.ArgumentParser()
parser.add_argument('chrom_sizes_path')
parser.add_argument('--cutwin', type=int, default=10)
args = parser.parse_args()
cutwin = args.cutwin

with open(args.chrom_sizes_path) as lines_in:
    chrom_sizes = dict([line.strip().split() for line in lines_in])
for key in chrom_sizes:
    chrom_sizes[key] = int(chrom_sizes[key]) - cutwin

with sys.stdin as lines_in:
    for line in lines_in:
        line = line.strip().split()
        refname = line[RNAME]
        tlen = int(line[TLEN])
        strand = '+' if int(line[FLAG]) & 16 else '-'
        if strand == '+':
            tn5_shift = 4
            cut_coord = (int(line[POS]) + tn5_shift) - 1
        else:
            tn5_shift = -5
            cut_coord = (int(line[POS]) + len(line[SEQ]) + tn5_shift) + 1
        cut_coord = max(min(cut_coord, chrom_sizes[refname]), 0)

        cut_out = '{!s}\t{!s}\t{!s}\t{!s}\t1000\t{!s}\n'.format(refname, max(cut_coord - cutwin, 0), min(cut_coord + cutwin, chrom_sizes[refname] + cutwin), line[QNAME]+'.pos_{!s}'.format(line[POS]), strand)
        sys.stdout.write(cut_out)
