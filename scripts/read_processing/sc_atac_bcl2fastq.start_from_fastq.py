#Python pipeline to convert NextSeq BCL files to fastq,
#fix errors in barcodes and trim adapters from fastqs.

#Current issues:
# - Hardcodes trimmomatic as located in Darren's
#   home directory.

import argparse
import os
from plumbum import local
import subprocess
import glob

TRIMM = 'src/Trimmomatic-0.36/trimmomatic-0.36.jar'
SRCDIR = os.path.dirname(__file__)

parser = argparse.ArgumentParser(description='A program to make fastq files from bcl2fastq ready for scATAC-seq analysis.')
parser.add_argument('--read1', help='Path to read1 fastq file.')
parser.add_argument('--read2', help='Path to read2 fastq file.')
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
#parser.add_argument('-S','--sample-sheet',help='Sample sheet location',dest='samplesheet',required=False,default="/net/Sample_sheet.csv")
parser.add_argument('-X','--nextseq',help='NextSeq run indicator',dest='nextseq',action="store_true")
#parser.add_argument('-H','--hiseq',help='HiSeq run indicator',dest='hiseq',action="store_true")
args = parser.parse_args()

##def submitter(commander):
##        """Submits commands directly to the command line and waits for the process to finish."""
##        submiting = subprocess.Popen(commander,shell=True)
##        submiting.wait()

#try:
#	os.makedirs(args.outdir)
#except OSError:
#	print 'Outdir already exists...'

#if not args.hiseq:
#	print "Making fastq files..."
#	print "Converting BCL files..."
#        local['bcl2fastq']('-R', args.rundir, '-o', os.path.join(args.outdir, 'fastq'), '--no-lane-splitting', '--sample-sheet', args.samplesheet)
##	fastqer = 'bcl2fastq -R ' + args.rundir + ' -o ' + os.path.join(args.outdir, 'fastq/') + ' --no-lane-splitting --sample-sheet ' + args.sampsheet
##	submitter(fastqer)

print "Fixing barcodes..."

splitter_params = [os.path.join(SRCDIR, 'NCP_fastq_10bpbarcode_split.py'), 
                   '-1', args.read1,
                   '-2', args.read2,
                   '-O1', os.path.join(args.outdir, args.prefix + '.split.1.fq'),
                   '-O2', os.path.join(args.outdir, args.prefix + '.split.2.fq'),
                   '-L', os.path.join(args.outdir, args.prefix + '.split.log'),
                   '-Z']
if args.nextseq:
        splitter_params.append('-X')
#splitter = 'python ' + os.path.join(SRCDIR, 'NCP_fastq_10bpbarcode_split.py') + ' -1 ' + os.path.join(args.outdir, 'fastq/Undetermined_S0_R1_001.fastq.gz') + ' -2 ' + os.path.join(args.outdir, 'fastq/Undetermined_S0_R2_001.fastq.gz') + ' -O1 ' + os.path.join(args.outdir, args.prefix + '.split.1.fq') + ' -O2 ' + os.path.join(args.outdir, args.prefix + '.split.2.fq'), ' -L ' + os.path.join(args.outdir, args.prefix + '.split.log') + ' -Z'
#if args.nextseq:
#	splitter = splitter + ' -X'

#if args.hiseq:
#        splitter_params = [os.path.join(SRCDIR, 'NCP_fastq_HiSeq_10bpbarcode_split.py'), 
#                           '-1', os.path.join(args.rundir, 'Undetermined_S0_R1_001.fastq.gz'),
#                           '-2', os.path.join(args.rundir, 'Undetermined_S0_R2_001.fastq.gz'),
#                           '-I1', os.path.join(args.rundir, 'Undetermined_S0_I1_001.fastq.gz'),
#                           '-I2', os.path.join(args.rundir, 'Undetermined_S0_I2_001.fastq.gz'),
#                           '-I3', os.path.join(args.rundir, 'Undetermined_S0_I3_001.fastq.gz'),
#                           '-I4', os.path.join(args.rundir, 'Undetermined_S0_I4_001.fastq.gz'),
#                           '-O1', os.path.join(args.outdir, args.prefix + '.split.1.fq'),
#                           '-O2', os.path.join(args.outdir, args.prefix + '.split.2.fq'),
#                           '-L' + os.path.join(args.outdir, args.prefix + '.split.log'),
#                           '-Z']
#	splitter = 'python ' + os.path.join(SRCDIR, 'NCP_fastq_HiSeq_10bpbarcode_split.py') + ' -1 ' + os.path.join(args.rundir, 'Undetermined_S0_R1_001.fastq.gz') + ' -2 ' + os.path.join(args.rundir, 'Undetermined_S0_R2_001.fastq.gz') + ' -I1 ' + os.path.join(args.rundir, 'Undetermined_S0_I1_001.fastq.gz') + ' -I2 ' + os.path.join(args.rundir, 'Undetermined_S0_I2_001.fastq.gz') + ' -I3 ' + os.path.join(args.rundir, 'Undetermined_S0_I3_001.fastq.gz') + ' -I4 ' + os.path.join(args.rundir, 'Undetermined_S0_I4_001.fastq.gz') + ' -O1 ' + os.path.join(args.outdir, args.prefix + '.split.1.fq') + ' -O2 ' + os.path.join(args.outdir, args.prefix + '.split.2.fq') + ' -L ' + os.path.join(args.outdir, args.prefix + '.split.log') + ' -Z'
#submitter(splitter)
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
#trimmer = 'java -Xmx1G -jar ' + TRIMM + ' PE ' + os.path.join(args.outdir, args.prefix + '.split.1.fq.gz') + ' ' + os.path.join(args.outdir, args.prefix + '.split.2.fq.gz') + ' ' + os.path.join(args.outdir, args.prefix + '.split.1.trimmed.paired.fastq.gz') + ' ' + os.path.join(args.outdir, args.prefix + '.split.1.trimmed.unpaired.fastq.gz') + ' ' + os.path.join(args.outdir, args.prefix + '.split.2.trimmed.paired.fastq.gz') + ' ' + os.path.join(args.outdir, args.prefix + '.split.2.trimmed.unpaired.fastq.gz') + ' ILLUMINACLIP:' + os.path.join(os.path.dirname(TRIMM), 'adapters/NexteraPE-PE.fa') + ':2:30:10:1:true TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 2> ' + os.path.join(args.outdir, args.prefix + '.split.trimmomatic.log')
#submitter(trimmer)

print "Cleaning up..."

local['rm'](os.path.join(args.outdir, args.prefix + '.split.1.fq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.2.fq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.1.trimmed.unpaired.fastq.gz'),
            os.path.join(args.outdir, args.prefix + '.split.2.trimmed.unpaired.fastq.gz'))
#local['gzip'](os.path.join(args.outdir, 'fastq/Undetermined*.fastq'))

#cleaner = 'rm ' + os.path.join(args.outdir, args.prefix + '.split.1.fq.gz') + '; rm ' + os.path.join(args.outdir, args.prefix + '.split.2.fq.gz') + '; rm ' + os.path.join(args.outdir, args.prefix + '.split.1.trimmed.unpaired.fastq.gz') + '; rm ' + os.path.join(args.outdir, args.prefix + '.split.2.trimmed.unpaired.fastq.gz') + '; gzip ' + os.path.join(args.outdir, 'fastq/Undetermined*.fastq')
#submitter(cleaner)
