#!/bin/sh

#1 - Reference genome path
reference=$1

# Reference genome folder
ref_folder=$(dirname $reference)

# Move to reference genome folder
cd $ref_folder

# Create index of reference genome
module load Bioinformatics
module load bwa

bwa index $reference

# Create fai index file
module purge
module load Bioinformatics 
module load samtools

samtools faidx $reference

# Create sequence dictionary
module purge
module load Bioinformatics
module load gatk

gatk CreateSequenceDictionary -R $reference