#!/bin/bash

# This is a pipeline for mapping exons onto transcripts and taking intersection of overlapping exons

n_threads=6
best_alignment=True
e_val=1e-10
word_size=9
match_score=1
mismatch_score=-1
gapopen=2
gapextend=1


echo 'FORMATTING DATABASE'
# Format database, i.e. file with gene models generated with 'make_gene_model_from_exons.py'

formatdb -pF -i ../ensembl_gene_model/gene_model.fasta

echo 'MAPPING EXONS TO REFERENCE'
# Map exons onto transcripts using distant megablast parameters

python map_seq_on_exons_by_blast.py -r ../ensembl_gene_model/gene_model.fasta -g ../ensembl_gene_model/gene_model.gff3 -q ../sample/transcripts_sample.fasta -b ../ensembl_gene_model/gene_model.fasta -f output_file -e $e_val -a $best_alignment -w $word_size -n $n_threads -r $match_score -s $mismatch_score -y $gapopen -x $gapextend

# The output file 'exons_alignment_by_blast_out_global.gff3' contains information about coordinates of model exons on transcripts.

echo 'TAKING INTERSECTION OF THE OVERLAPING EXONS'
# Take intersection of overlapping exons (e.g. two overlapping exons will be split into three non-overlapping regions)

python exon_mip_intersection.py exons_alignment_by_blast_out_global.gff3

# The output file, exon_mip_intersection.out is a table containing three columns: reference transcript name, exon start position, exon stop position
