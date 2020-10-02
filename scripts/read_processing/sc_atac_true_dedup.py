#This script was originally written by
#Darren Cusanovich. It has been slightly modified
#to ensure that deduplication preserves mate pairs.

import pysam #v0.8.1
import sys

inbam = sys.argv[1]
outbam = sys.argv[2]

readsin = pysam.AlignmentFile(inbam,"rb")
readsout = pysam.AlignmentFile(outbam,"wb",template=readsin)
refs = readsin.references
#I just break it up by chromosome to give a status update AND to make sure the 'readdic' doesn't get too big, which might hinder speed
for refchrom in refs:
        if 'chrM' in refchrom or 'chrGL' in refchrom or 'chrNC' in refchrom or 'chrhs' in refchrom or 'random' in refchrom or 'chrU' in refchrom:
                continue
        readdic = {}
        print("Deduplicating " + refchrom + "...")
        for read in readsin.fetch(refchrom):
                readname = read.qname.split(':')[0]
                if 'CTF' in readname or 'AMBIG' in readname:
                        continue

                if read.tlen < 0:
                        fragstart = read.mpos
                        fragend = read.mpos - read.tlen #read.tlen is negative for reverse strand reads
                else:
                        fragstart = read.pos
                        fragend = read.pos + read.tlen
                keytup = (readname, fragstart, fragend)

                try:
                        chosen_mate = readdic[keytup]
                except KeyError:
                        readdic[keytup] = read
                        readsout.write(read)
                else:
                        if chosen_mate.qname == read.qname:
                                readsout.write(read)
