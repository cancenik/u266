#!/bin/bash -l

#SBATCH -t 24:00:00
#SBATCH --mem=4000
#SBATCH --account=mpsnyder

source $HOME/.bash_profile

## CC Update Jul 5 ##
## Changed variable handling to make the script compatible with SLURM


################# ADAPTER REMOVAL
# Remove 3' adapter using cutadapt
# -g 5'adapter -a 3'adapter -e ERROR_RATE (0.1)

# Clontech suggested adaptor
# Keeping the As to a stretch of 10 is a good idea as there are reads that clearly end with As but sequenced incorrectly after a longstretch of As. 
#cutadapt -m 15 -u 3 -a AAAAAAAAAA input.fastq > output.fastq
# Min length can be taken 15 as well. 
echo "cutadapt -a $2 --overlap=4 --trimmed-only --minimum-length=17 --quality-cutoff=33 $1 >Clipped_$1 2>>LOG_$1" >>LOG_$1
cutadapt -u 3 -a $2 --overlap=4 --trimmed-only --maximum-length=44 --minimum-length=15 --quality-cutoff=33 $1 >Clipped_$1 2>>LOG_$1 

# Example
# EMI ADAPTER ATCTCGTATGCCGTCTTCTGCTTG
# cutadapt -a CTGTAGGCACCATCAAT  --overlap=5 --trimmed-only --minimum-length=18 --quality-cutoff=33 
#################
################# ALIGNMENT
echo "bowtie2 -L 15 -x rRNA -q Clipped_$1  --un Unaligned_rRNA_$1 2>>LOG_$1 | samtools view -bS - >rRNA_aln.bam" >>LOG_$1
bowtie2 -L 15 -x rRNA -q Clipped_$1  --un Unaligned_rRNA_$1 2>>LOG_$1 | samtools view -bS - >rRNA_aln.bam
echo "bowtie2 -L 15 --norc -x APPRIS -q Unaligned_rRNA_$1 -S appris_aln.sam --un Unaligned_Appris_$1 2>>LOG_$1" >>LOG_$1
bowtie2 -L 15 --norc -x APPRIS -q Unaligned_rRNA_$1  --un Unaligned_Appris_$1 2>>LOG_$1 | samtools view -bS - >appris_aln.bam
#################

############# READ COUNTS & ANALYSIS 
## Convert SAM OUTPUT TO BAM FILE CUTOFF WITH -q
## OR Keep the results in sam format
samtools view -b -q 2 appris_aln.bam >appris_aln_Q2.bam

## Count reads per feature using BedTools force strandedness
## Need to specify --split when junction mapping is possible
coverageBed -b appris_aln_Q2.bam -a /srv/gsfs0/projects/snyder/ccenik/SEQUENCE_INDEX/RAW_SEQUENCES/GENCODE/Appris_Regions.bed -s -split >Transcript_Counts_$1

## OR MANUALLY COUNT FOR TRANSCRIPTS, Remove Reads with a non-zero samtools flag
#samtranscriptaln_tocounts.pl --file appris_aln_Q2.sam



## FASTQC AFTER ADAPTER REMOVAL
#echo "fastqc Clipped_$1" >>LOG_$1
#fastqc Clipped_$1
## FASTQC AFTER RRNA ALIGNMENT
#echo "fastqc Unaligned_rRNA_$1" >>LOG_$1
#fastqc Unaligned_rRNA_$1
## FASTQC TRANS ALIGNMENTS
#echo "appris_aln_Q2.bam" >>LOG_$1
#fastqc appris_aln_Q2.bam



