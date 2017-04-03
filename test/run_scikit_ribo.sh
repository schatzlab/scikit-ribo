#!/bin/sh
set -e

# --------------------
# Example shell script
# Run scikit-ribo
# --------------------

gtf=/gtf-file-location/
fasta=/ref-genome-fasta-location/
bam=/bam-file-location/ # from STAR
rnafold=/folder-to-rnafold-results/ # folder that contains *dp_ps files
prefix=/index-prefix-to-use/
index_folder=/index-output-folder/
RNA=/RNAseq-TPM-file-location/ # salmon or kallisto
output=/output-folder/
unmap=/unmap-regions/

scikit-ribo-build.py \
-g $gtf \
-f $fasta \
-p $prefix \
-r $rnafold \
-t $RNA \
-o $index_folder

scikit-ribo-run.py
-i $bam \
-f $index_folder \
-p $prefix \
-o $output \
-u $unmap
