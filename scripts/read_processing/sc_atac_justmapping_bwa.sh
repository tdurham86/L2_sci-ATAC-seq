#! /usr/bin/env bash

#This script is was originally
#written by Darren Cusanovich, and has
#been lightly modified to align reads to
#the ce10 and ce11 genomes with bowtie2.

if [ $# -ne 2 ]
then
        echo -e \\nWarning!!! Incorrect number of parameters!!!\\n\\nUsage: bash sc_atac_justmapping_bwa.sh [Work Prefix] [Species \(e.g. dm3 or hg19\)];
        echo -e;
        exit;
fi

if [ $2 = "hg19" ]
then
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol7/cusanovi/genomes/$2/hs37d5.fa $1.split.1.trimmed.paired.fastq.gz > $1.split.1.sai&
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol7/cusanovi/genomes/$2/hs37d5.fa $1.split.2.trimmed.paired.fastq.gz > $1.split.2.sai&
	#wait
	#~/bin/bwa/bwa sampe /net/shendure/vol7/cusanovi/genomes/$2/hs37d5.fa $1.split.1.sai $1.split.2.sai $1.split.1.trimmed.paired.fastq.gz $1.split.2.trimmed.paired.fastq.gz | samtools view -Sb - > $1.split.bam;
    bowtie2 -p $NSLOTS -X 2000 -3 1 -x proj/hnd-1_rnaseq/ref_data/hg19_just_base_chrs/hg19 -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam;
    samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrGL | grep -v '[0-9]'$'\t'chrNC | grep -v '[0-9]'$'\t'chrhs | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -T ./ -o $1.split.q10.sort.bam -@ $NSLOTS -;
fi

if [ $2 = "mm9" ]
then
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol7/cusanovi/genomes/$2/$2.fa $1.split.1.trimmed.paired.fastq.gz > $1.split.1.sai&
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol7/cusanovi/genomes/$2/$2.fa $1.split.2.trimmed.paired.fastq.gz > $1.split.2.sai&
	#wait
	#~/bin/bwa/bwa sampe /net/shendure/vol7/cusanovi/genomes/$2/$2.fa $1.split.1.sai $1.split.2.sai $1.split.1.trimmed.paired.fastq.gz $1.split.2.trimmed.paired.fastq.gz | samtools view -Sb - > $1.split.bam;
	bowtie2 -p 16 -X 2000 -3 1 -x /net/shendure/vol10/projects/scATAC/nobackup/genomes/mm9/tophat/mm9 -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam;
	samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chr'*'_random$'\t' | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 16 - $1.split.q10.sort;
fi

if [ $2 = "dm3" ]
then
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol10/projects/scATAC/nobackup/genomes/$2/$2.fa $1.split.1.trimmed.paired.fastq.gz > $1.split.1.sai&
	#~/bin/bwa/bwa aln -t 12 /net/shendure/vol10/projects/scATAC/nobackup/genomes/$2/$2.fa $1.split.2.trimmed.paired.fastq.gz > $1.split.2.sai&
	#wait
	#~/bin/bwa/bwa sampe /net/shendure/vol10/projects/scATAC/nobackup/genomes/$2/$2.fa $1.split.1.sai $1.split.2.sai $1.split.1.trimmed.paired.fastq.gz $1.split.2.trimmed.paired.fastq.gz | samtools view -Sb - > $1.split.bam;
	bowtie2 -p 16 -X 2000 -3 1 -x /net/shendure/vol10/projects/scATAC/nobackup/genomes/dm3/dm3 -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam;
	samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrU | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 16 - $1.split.q10.sort;
fi

if [ $2 = "ce10" ]
then
    bowtie2 -p $NSLOTS -X 2000 -3 1 -x ATAC_sequencing/2018_worm_atac/ref_data/WS230/c_elegans.WS230.genomic -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam;
    samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v 'MtDNA' | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -T ./ -o $1.split.q10.sort.bam -@ $NSLOTS -;
fi

if [ $2 = "ce11" ]
then
    bowtie2 -p $NSLOTS -X 2000 -3 1 -x ATAC_sequencing/2018_worm_atac/ref_data/WS235/c_elegans.WS235.genomic -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam;
    samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v 'MtDNA' | grep -v 'chrM' | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -T ./ -o $1.split.q10.sort.bam -@ $NSLOTS -;
fi

samtools index $1.split.q10.sort.bam;

python ATAC_sequencing/all_sciATAC/src/preprocess/sc_atac_true_dedup.py $1.split.q10.sort.bam $1.true.nodups.bam;

samtools view $1.true.nodups.bam | sort -u -k1,1 | cut -f9 > $1.insertsize.txt;

samtools index $1.true.nodups.bam;
