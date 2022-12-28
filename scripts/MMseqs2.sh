#!/bin/bash

OUT=$1
CDS=$2
ID=$3

###################### concatenate input CDS files ########################

cat $CDS/*.fna > $OUT/concatenated_CDSes.fna

######################### Create query and target databases ########################

echo "Creating query and target databases from input CDS files..."

# Query database
mmseqs createdb $OUT/concatenated_CDSes.fna $OUT/QueryDB

# Target database
mmseqs createdb $OUT/concatenated_CDSes.fna $OUT/TargetDB

######################### Create folder for temporary files ########################
mkdir $OUT/tmp

######################### Index the target database ################################
mmseqs createindex $OUT/TargetDB $OUT/tmp --search-type 3

######################### Search for sequence homology #############################

echo "Searching for sequence homology using MMseqs2..."

mmseqs search -a --min-seq-id $ID $OUT/QueryDB $OUT/TargetDB $OUT/ResultDB $OUT/tmp

######################### Convert output file #######################################

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