from optparse import OptionParser
import os
import time
import sys
from Bio import SeqIO

formats = ['fastq', 'fasta', 'fa', 'fas']


class Exon:
    """ Data storage for Exons. """

    def __init__(self, name, o_start, o_stop, strand, seq):
        self.name = name
        self.o_start = o_start
        self.o_stop = o_stop
        self.strand = strand
        self.seq = seq
        self.length = len(seq)
        self.start_five_p = 0
        self.stop_three_p = 0

        self.place_in_gene = False


class Gene:
    """ Data storage and processor for genes. """

    def __init__(self, name, chrom_name):

        self.name = name
        self.chrom_name = chrom_name
        self.exons_dic = {}
        self.strand = 1
        self.gene_template_list = []
        self.gene_template_dic = {}
        self.exons_list = []
        self.template_length = 0

        self.min_coord = 0
        self.max_coord = 0

    def add_exon(self, exon_name, start, stop, strand, seq):

        ex = Exon(exon_name, int(start), int(stop), int(strand), seq)

        if exon_name not in self.exons_dic:
            self.exons_dic[exon_name] = ex
            self.strand = int(strand)

    def get_exons_list(self):
        return self.exons_dic.values()

    def gene_template_with_strand_corrections(self):
        """ Generate gene template; adjust data for genes on Reverse strand. """

        self.exons_list = self.get_exons_list()

        # change of coordinate system by mirroring coordinates with reflection plane on side of the highest coord.
        if self.strand == -1:

            the_highest_coordinate_of_all_exons = max([ex.o_stop for ex in self.exons_list])
            the_highest_coordinate_of_all_exons_shifted = the_highest_coordinate_of_all_exons + 100

            # change exons coords:
            for _, ex in self.exons_dic.iteritems():
                start = ex.o_start
                stop = ex.o_stop
                ex.o_start = stop + 2 * (the_highest_coordinate_of_all_exons_shifted - stop)
                ex.o_stop = start + 2 * (the_highest_coordinate_of_all_exons_shifted - start)

        self.exons_list = self.exons_dic.values()

        ex_tup_list = []
        new_exons_list = []

        for ex in self.exons_list:
            ex_tup_list.append((ex.name, ex.o_start))

        local_sorted_exons_list = sorted(ex_tup_list, key=lambda length: length[1], reverse=False)

        for ex_in_tuple in local_sorted_exons_list:
            new_exons_list.append(self.exons_dic[ex_in_tuple[0]])

        # shifts all coordinates to 1 for first base
        first_coord = new_exons_list[0].o_start
        for ex in new_exons_list:
            ex.start_five_p = ex.o_start - first_coord + 1
            ex.stop_three_p = ex.o_stop - first_coord + 1

        max_min = []
        for ex in self.exons_list:
            max_min.append((ex.o_start, ex.o_stop))

        max_min_list = sorted(max_min, key=lambda length: length[0], reverse=False)
        self.min_coord = max_min_list[0][0]

        max_min_list = sorted(max_min, key=lambda length: length[1], reverse=True)
        self.max_coord = max_min_list[0][1]

        self.exons_list = new_exons_list

        self.gene_template_dic = {}
        self.gene_template_list = []

        # coordinates must start from 1
        for i in xrange(self.max_coord - self.min_coord + 2):
            self.gene_template_dic[i] = "-"
            self.gene_template_list.append("-")

        self.template_length = len(self.gene_template_list)

    def generate_gene_model_seq(self):
        """ Generates gene model sequence. """

        if len(self.gene_template_list) == 0:
            raise Exception("Problem in gene_template_list")

        return str(''.join([ex[0] for ex in self.gene_template_list]))

    def generate_gene_info_for_exons(self):
        """ Generates dictionary output with exons annotation in gene. """

        if len(self.gene_template_list) == 0:
            raise Exception("Problem in gene_template_list")

        exons_dic = {}

        for idx in xrange(len(self.gene_template_list)):

            seq_idx = idx + 1
            exon_name = self.gene_template_list[idx][1]
            if exon_name not in exons_dic:
                exons_dic[exon_name] = [seq_idx, seq_idx]
            else:
                start, stop = exons_dic[exon_name]
                exons_dic[exon_name] = [start, seq_idx]

        return exons_dic

    def generate_model(self):
        """ Generates gene model - sequence & annotation info. """

        gene_template_list_l = ['-' for i in xrange(self.max_coord + 2)]

        for i in xrange(len(self.exons_list)):
            exon = self.exons_list[i]

            idx_in_seq = 0
            for base_idx in range(exon.o_start, exon.o_stop + 1):
                gene_template_list_l[base_idx] = (exon.seq[idx_in_seq], exon.name)
                idx_in_seq += 1

        self.gene_template_list = [x for x in gene_template_list_l if x != '-']

        exons_list_to_gff = []
        exons_dic = self.generate_gene_info_for_exons()
        for exon in self.exons_list:
            if exon.name not in exons_dic: continue

            start, stop = exons_dic[exon.name]
            exons_list_to_gff.append([exon.name, start, stop])

        return self.generate_gene_model_seq(), exons_list_to_gff


def check_format(name):
    """ Simple checker for sequence input file.  """

    for f in formats:
        if str(name).endswith(f):
            return True
    return False


def main(args=[]):
    """ Main function to process data in batch mode, with description of each stage. """

    # --- Argument parsing

    usage = "usage: %prog [options] arg \nProgram generate gene models for exons' dataset from ENSEMBL"
    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-f", "--fasta_file", dest="FASTA_FILE", help="fasta file", action="store", type="string")
    parser.add_option("-e", "--exons_info", dest="EXONS", help="exons description information", action="store",
                      type="string")
    parser.add_option("-o", "--output_folder", dest="OUTPUT_FOLDER", help="output folder")

    (options, arg) = parser.parse_args(args)

    # ---  Entering program

    sys.stdout.write("\nEntering program\n")
    t_st = time.time()

    if not os.path.isdir(options.OUTPUT_FOLDER):
        raise Exception('\nWrong output directory!')

    if not check_format(options.FASTA_FILE):
        raise Exception('\nWrong input fasta file! This is not fasta format')

    os.chdir(options.OUTPUT_FOLDER)
    log_info_hlr = open(options.OUTPUT_FOLDER + os.sep + "outoutinfo.log", "w")

    # ---  parsing fasta file

    exons_sequences_dic = {}

    log_str = "Parsing: {}\n".format(options.FASTA_FILE)
    log_info_hlr.write(log_str)
    sys.stdout.write(log_str)

    for seq in SeqIO.parse(open(options.FASTA_FILE, 'r'), "fasta"):
        exons_sequences_dic[seq.id] = str(seq.seq)

    # ---   parsing exons' info file

    # exons' info schema:
    # Chromosome Name    Ensembl Gene ID    Strand    Ensembl Exon ID    Exon Chr Start (bp)    Exon Chr End (bp)
    # GL174807.1    ENSXETG00000025632    -1    ENSXETE00000426649    14618    15080
    # GL174807.1    ENSXETG00000025632    -1    ENSXETE00000397496    13786    13988

    log_str = "Parsing: {}\n".format(options.EXONS)
    log_info_hlr.write(log_str)
    sys.stdout.write(log_str)

    genes_dic = {}
    desc_exons_hlr = open(options.EXONS, 'r')

    for line in desc_exons_hlr:

        if len(str(line).strip()) == 0 or str(line).strip().startswith("#"): continue

        chr_name = str(line).strip().split('\t')[0]
        gene = str(line).strip().split('\t')[1]
        strand = str(line).strip().split('\t')[2]
        exon = str(line).strip().split('\t')[3]
        ex_start = str(line).strip().split('\t')[4]
        ex_end = str(line).strip().split('\t')[5]

        if exon not in exons_sequences_dic:
            log_info_hlr.write("\nNo exon sequence: " + str(exon))
            continue

        exon_seq = exons_sequences_dic[exon]

        if gene not in genes_dic:
            g = Gene(gene, chr_name)
            genes_dic[gene] = g

        g = genes_dic[gene]
        g.add_exon(exon, ex_start, ex_end, strand, exon_seq)

    desc_exons_hlr.close()

    # ---   generating gene models (fasta and coordinates/info file)

    gene_model_seq_hlr = open(options.OUTPUT_FOLDER + os.sep + "gene_model" + ".fasta", 'w')
    gene_model_gff3_hlr = open(options.OUTPUT_FOLDER + os.sep + "gene_model" ".gff3", 'w')

    gene_model_counter = 0
    for gene_name, gene_obj in genes_dic.iteritems():

        gene_model_counter += 1
        log_str = "Generating gene model ({}): {}\n".format(gene_model_counter, gene_name)
        sys.stdout.write(log_str)

        # correct gene object
        gene_obj.gene_template_with_strand_corrections()

        # make gene model!
        seq_string, exons_list_to_gff = gene_obj.generate_model()

        gene_model_seq_hlr.write(">" + gene_name + "\n")
        gene_model_seq_hlr.write(seq_string + "\n")

        for exon in exons_list_to_gff:
            group = "Name=" + str(exon[0]) + ";" + "Last=" + str(exon == exons_list_to_gff[len(exons_list_to_gff) - 1])
            string_gff3 = gene_obj.name + "\t" + "." + "\t" + "gene_model_exons" + "\t" + str(
                exon[1]) + "\t" + str(exon[2]) + "\t" + "." + "\t" + str(
                gene_obj.strand) + "\t" + "." + "\t" + group
            gene_model_gff3_hlr.write(string_gff3 + "\n")

        genes_dic[gene_name] = None

    # ---   closing program

    gene_model_seq_hlr.close()
    gene_model_gff3_hlr.close()

    t_end = time.time()

    sys.stdout.write("\nwork done...")
    sys.stdout.write("\nprocess time [s]: " + str(t_end - t_st))

    log_info_hlr.close()


if __name__ == "__main__":
    main(sys.argv[1:])
