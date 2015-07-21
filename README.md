# targeted-resequencing-with-mips


1. Download from ensembl/biomart exon sequences and information about exon coordinates of Xenopus, humans, chicken and anolis genes of interest
two files, e.g. 'xenopus69_immuno_exons.fasta' and 'xenopus69_immuno_exons_info.txt'

2. use 'make_gene_model_from_exons.py' to construct gene models (exons from one gene merged together) and gff3 file (coordinates of exons on these gene models)
files: 'gene_model.fasta' & 'gene_model.gff3'
merge all fasta and gff3 files of all species together
(can check with 'exon_stats.py' how many overlapping exons there are)


now run script: find_exons.sh
3. filter gff3 file to remove redundant exons which share at least one identical coordinate with 'filter_exons.py'
(can check with 'exon_stats.py' how many overlapping exons there are after filtering)

4. use the gff3 file from the output as a gene model file
index database and run 'map_seq_on_exons_by_blast.py'
parameters of the blastn search are (distant search):

```
e_val=1e-6
word_size=9
match_score=1
mismatch_score=-1
gapopen=2
gapextend=1
```

5. the resulting file 'exons_alignment_by_blast_out_global.gff3' has coordinates of model exons put onto reference contigs,
check with exon_stats.py how many overlapping exons there are (careful with the last column and location of exonID)

6. use 'exon_intersection.py' to filter gff3 file where different regions hit by the same exon are merged and intersection is taken from the overlapping exons
output 'exon_intersection.out' is a table, with ref name in the first column, real start and stop positions in the second and third, and some info in the fourth column




