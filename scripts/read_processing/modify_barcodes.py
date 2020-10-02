#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
from plumbum import local
import sys

parser = argparse.ArgumentParser()
parser.add_argument('bam_path')
parser.add_argument('append_to_barcode')
args = parser.parse_args()

bam_lines = local['samtools']['view', '-h', args.bam_path]
for line in bam_lines.popen().iter_lines():
    if line[0][0] != '@':
        bc = line[0].rstrip().split(':')
        bc[0] = bc[0] + args.append_to_barcode
        sys.stdout.write(':'.join(bc) + '\n')
    else:
        sys.stdout.write(line[0].rstrip() + '\n')
