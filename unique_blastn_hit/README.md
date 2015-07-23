
### This step is to select transcripts which unambiguously map only to one model species gene

* First prepare a file of blast results in 24 col NCBI format

    #### Example
    
    Here is an example distant megablast with same parameters that were used for mapping exons to transcripts
    
    ```
    blastn -db ../ensembl_gene_model/gene_model.fasta -query ../sample/transcripts_sample.fasta -evalue 1e-10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen" -out blast_output -word_size 9 -num_threads 6 -reward 1 -penalty -1 -gapopen 2 -gapextend 1 -num_alignments 10

    ```
* Run *Process_blast_output_24_col.py* to select only unique pairs transcript-gene model

    ```
    python Process_blast_output_24_col.py blast_output
    ```
    
    Unique hits are in *Single_hit.txt* file