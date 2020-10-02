#This script is a lightly modified version of code 
#originally written by Darren Cusanovich. It has 
#been adapted to use the plumbum package for
#executing shell commands.

#Python pipeline to convert NextSeq BCL files to fastq,
#fix errors in barcodes and trim adapters from fastqs.

#Current issues:
# - Hardcodes trimmomatic jar location.

import argparse
import os
from plumbum import local
import subprocess
import glob

TRIMM = 'src/Trimmomatic-0.36/trimmomatic-0.36.jar'
SRCDIR = os.path.dirname(__file__)

parser = argparse.ArgumentParser(description='A program to convert NextSeq BCL files to fastq files for scATAC-seq analysis.')
parser.add_argument('-R','--rundir', help='Run directory containing BCL files',dest='rundir',required=True)
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
parser.add_argument('-S','--sample-sheet',help='Sample sheet location',dest='samplesheet',required=False,default="/net/Sample_sheet.csv")
parser.add_argument('-X','--nextseq',help='NextSeq run indicator',dest='nextseq',action="store_true")
parser.add_argument('-H','--hiseq',help='HiSeq run indicator',dest='hiseq',action="store_true")
args = parser.parse_args()

try:
	os.makedirs(args.outdir)
except OSError:
	print 'Outdir already exists...'

if not args.hiseq:
	print "Making fastq files..."
	print "Converting BCL files..."
        local['bcl2fastq']('-R', args.rundir, '-o', os.path.join(args.outdir, 'fastq'), '--no-lane-splitting', '--sample-sheet', args.samplesheet)

print "Fixing barcodes..."

splitter_params = [os.path.join(SRCDIR, 'NCP_fastq_10bpbarcode_split.py'), 
                   '-1', os.path.join(args.outdir, 'fastq/Undetermined_S0_R1_001.fastq.gz'),
                   '-2', os.path.join(args.outdir, 'fastq/Undetermined_S0_R2_001.fastq.gz'),
                   '-O1', os.path.join(args.outdir, args.prefix + '.split.1.fq'),
                   '-O2', os.path.join(args.outdir, args.prefix + '.split.2.fq'),
                   '-L', os.path.join(args.outdir, args.prefix + '.split.log'),
                   '-Z']
if args.nextseq:
        splitter_params.append('-X')

if args.hiseq:
        splitter_params = [os.path.join(SRCDIR, 'NCP_fastq_HiSeq_10bpbarcode_split.py'), 
                           '-1', os.path.join(args.rundir, 'Undetermined_S0_R1_001.fastq.gz'),
                           '-2', os.path.join(args.rundir, 'Undetermined_S0_R2_001.fastq.gz'),
                           '-I1', os.path.join(args.rundir, 'Undetermined_S0_I1_001.fastq.gz'),
                           '-I2', os.path.join(args.rundir, 'Undetermined_S0_I2_001.fastq.gz'),
                           '-I3', os.path.join(args.rundir, 'Undetermined_S0_I3_001.fastq.gz'),
                           '-I4', os.path.join(args.rundir, 'Undetermined_S0_I4_001.fastq.gz'),
                           '-O1', os.path.join(args.outdir, args.prefix + '.split.1.fq'),
                           '-O2', os.path.join(args.outdir, args.prefix + '.split.2.fq'),
                           '-L' + os.path.join(args.outdir, args.prefix + '.split.log'),
                           '-Z']

local['python'].__getitem__(splitter_params)()

print "Trimming adapters..."

#run Trimmomatic
with open(os.path.join(args.outdir, args.prefix + '.split.trimmomatic.log'), 'w') as out:
          local['java']('-Xmx1G', '-jar', TRIMM, 'PE', 
                        os.path.join(args.outdir, args.prefix + '.split.1.fq.gz'),
                        os.path.join(args.outdir, args.prefix + '.split.2.fq.gz'),
                        os.path.join(args.outdir, args.prefix + '.split.1.trimmed.paired.fastq.gz'),
                        os.path.join(args.outdir, args.prefix + '.split.1.trimmed.unpaired.fastq.gz'),
                        os.path.join(args.outdir, args.prefix + '.split.2.trimmed.paired.fastq.gz'),
                        os.path.join(args.outdir, args.prefix + '.split.2.trimmed.unpaired.fastq.gz'),
                        'ILLUMINACLIP:{!s}:2:30:10:1:true'.format(os.path.join(os.path.dirname(TRIMM), 'adapters/NexteraPE-PE.fa')), 
                        'TRAILING:3', 'SLIDINGWINDOW:4:10', 'MINLEN:20',
                        stderr=out)

print "Cleaning up..."

local['rm'](os.path.join(args.outdir, args.prefix + '.split.1.fq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.2.fq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.1.trimmed.unpaired.fastq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.2.trimmed.unpaired.fastq.gz'))
