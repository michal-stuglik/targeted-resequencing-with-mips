# targeted-resequencing-with-mips

#### Pipeline for identifying exon boundaries in transcripts using pairwise alignments with gene models of model species

## Construction of gene models of model species

Available scripts are in *ensembl_gene_model* directory

* Download from ensembl/biomart exon sequences and information about exon coordinates of model species genes. Example of two files for *Xenopus sp.*  can be found in *sample* directory

* Run *make_gene_model_from_exons.py* to construct gene models (exons from one gene will be merged together) and gff3 file (coordinates of exons on these gene models).

    ```
    python make_gene_model_from_exons.py -f xenopus_exons.fasta -e xenopus_exons_info.txt -o ./
    ```
    
    There are two output files: *gene_model.fasta* with exon sequences of one gene merged together and *gene_model.gff3* with coordinates of exons on these gene models.



## Mapping transcripts to gene models

Available scripts are in *exon_boundaries_on_transcript* directory. In this part exons are mapped to transcripts and exon boundaries on transcripts are identified. Run all the steps using *find_exons.sh* bash script.

### Example

* Format database, i.e. file with gene models

    ```
    formatdb -pF -i gene_model.fasta
    ```

* Map exons onto transcripts using distant megablast parameters

    ```
    n_threads=6
    best_alignment=True
    e_val=1e-10
    word_size=9
    match_score=1
    mismatch_score=-1
    gapopen=2
    gapextend=1

    python map_seq_on_exons_by_blast.py -r gene_model.fasta -g gene_model.gff3 -q transcript_sequences.fasta -b gene_model.fasta -f output_file -e $e_val -a $best_alignment -w $word_size -n $n_threads -r $match_score -s $mismatch_score -y $gapopen -x $gapextend
    ```
    The output file *exons_alignment_by_blast_out_global.gff3* contains information about coordinates of model exons on transcripts.
    
* Check overlapping exons with *exon_stats.py*.
    
    ```
    python exon_stats.py exons_alignment_by_blast_out_global.gff3
    ```

* Take intersection of overlapping exons (e.g. two overlapping exons will be split into three non-overlapping regions)

    ```
    python exon_mip_intersection.py
    ```
    
## Selecting unique blastn transcript-gene model pairs

Available script can be found in *unique_blastn_hit* directory. This step is to select only transcripts which unambiguously map only to one model species gene. 

First prepare a file of blast results in 24 col NCBI format. Then run *process_blast_output_24_col.py*
```
python process_blast_output_24_col.py blast_output
```



