#!/bin/bash

gene_model_file=gene_model_new.gff3
reference_file=mapping_mips_references.fasta
n_threads=6
best_alignment=True
e_val=1e-10
word_size=9
match_score=1
mismatch_score=-1
gapopen=2
gapextend=1

output_file=exons
database=gene_model_new.fasta


echo 'FILTERING GENE MODEL GFF3 FILE'

python filter_exons.py $gene_model_file

mv filter_exons.out gene_model_filtered.gff3
rm filter_exons.out

echo 'FORMATTING DATABASE'

formatdb -pF -i $database

echo 'MAPPING EXONS TO REFERENCE'

python map_seq_on_exons_by_blast.py -r $database -g gene_model_filtered.gff3 -q $reference_file -b $database -f $output_file -e $e_val -a $best_alignment -w $word_size -n $n_threads -r $match_score -s $mismatch_score -y $gapopen -x $gapextend

echo 'STATISTICS FOR REFERENCE GFF3 FILE'

python exon_stats.py exons_alignment_by_blast_out_global.gff3

echo 'TAKING INTERSECTION OF THE OVERLAPING EXONS'

python exon_mip_intersection.py


