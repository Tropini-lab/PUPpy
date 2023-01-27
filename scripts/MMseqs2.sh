#!/bin/bash

# Script: MMseqs2.sh
# Part of PUPpy distribution.
# Run from 'puppy-align' script
# Runs MMseqs2 search workflow to perform pairwise local alignments between genes in the defined microbial community.

# Author = Hans Ghezzi
# Credits = Hans Ghezzi
# Licence = GPL3
# Version = 1.0.0
# Maintainer = Hans Ghezzi
# Email = hans.ghezzi@gmail.com
# Status = Development

# Import user-defined variables from puppy-align script. These are the values of the 3 flags of the puppy-align script.
OUT=$1 # outdir defined in puppy-align
CDS=$2 # CDSdire defined in puppy-align
ID=$3 # ID% threshold defined in puppy-align

# Concatenate input CDS files
cat $CDS/*.fna > $OUT/concatenated_CDSes.fna

# Create query and target databases
echo "Creating query and target databases from input CDS files..."

# Create query database from the concatenated CDSes
mmseqs createdb $OUT/concatenated_CDSes.fna $OUT/QueryDB

# Create target database from the concatenated CDSes
mmseqs createdb $OUT/concatenated_CDSes.fna $OUT/TargetDB

# Create folder for temporary files
mkdir $OUT/tmp

# Index the target database
mmseqs createindex $OUT/TargetDB $OUT/tmp --search-type 3

# Search for sequence homology between query and target database. Alignments with ID% below threshold are not reported.
echo "Searching for sequence homology using MMseqs2..."
mmseqs search -a --min-seq-id $ID $OUT/QueryDB $OUT/TargetDB $OUT/ResultDB $OUT/tmp

# Convert output file
echo "Converting alignments to output file..."
mmseqs convertalis $OUT/QueryDB $OUT/TargetDB $OUT/ResultDB $OUT/ResultDB.tsv --format-output "query,target,qlen,tlen,alnlen,qstart,qend,tstart,tend,pident,qcov,tcov,evalue" --search-type 3

# Clean output directory
mv $OUT/tmp $OUT/mmseqs_tmp
mkdir -p $OUT/tmp
mv $OUT/Query* $OUT/tmp
mv $OUT/Target* $OUT/tmp
mv $OUT/Result* $OUT/tmp
mv $OUT/*concat* $OUT/tmp
mv $OUT/tmp/ResultDB.tsv $OUT