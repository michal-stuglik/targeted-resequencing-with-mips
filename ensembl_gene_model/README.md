### Gene models based on ENSEMBL exon dataset

Python scripts to generate gene models in 5'->3' orientation, based on info for all exons extracted from Biomart ENSEMBL (see sample folder for example files)


### Example

* Download from Biomart ENSEMBL exon sequences and information about exon coordinates of model species genes. Example of two files for *Xenopus sp.*  can be found in *sample* directory.
    First file is a simple fasta file with exon sequences and Exon ID in the sequence header. The second file consists of six columns:
    * Chromosome Name
    * Ensembl Gene ID
    * Strand
    * Ensembl Exon ID
    * Exon Chr start (bp)
    * Exon Chr stop (bp)

* Run *make_gene_model_from_exons.py* to construct gene models (exons from one gene will be merged together) and gff3 file (coordinates of exons on these gene models).

    ```
    python make_gene_model_from_exons.py -f ../sample/xenopus_exons.fasta -e ../sample/xenopus_exons_info.txt -o ./
    ```
    
    There are two output files: *gene_model.fasta* with exon sequences of one gene merged together and *gene_model.gff3* with coordinates of exons on these gene models.


