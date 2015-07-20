from optparse import OptionParser
import os
import time
import sys
from Bio.Blast import NCBIXML

formats = ['fastq', 'fasta', 'fa', 'fas']


class AlignmentBase:
    """ Data structure to store alignment object. """

    def __init__(self, subject_base, query_base, position_subject, position_query, position_query_local, query_name, subject_name, score, strand, query_geneName_local):
        self.s = subject_base
        self.q = query_base

        self.position_subject = position_subject
        self.position_query = position_query
        self.position_query_local = position_query_local
        self.exons_list = None

        self.query_name = query_name
        self.subject_name = subject_name
        self.score = score
        self.strand = strand
        self.query_geneName_local = query_geneName_local

    def set_exons(self, exons_list):
        self.exons_list = exons_list


def save_exoninfo_in_gff(out_gff3_hlr, *args):
    out_gff3_hlr.write("\t".join(args) + "\n")


def extract_and_write_gff(alignments_objects_list, out_gff3_global_hlr, out_gff3_local_hlr):

    gff3_source = "."
    gff3_frame = "."
    gff3_feature = "blast_hit_model"

    exons_dic = {}  # exon_name: [position_query, position_query_local, last_exon_info]

    for al_obj_in_hsp in alignments_objects_list:
        for exon_info in al_obj_in_hsp.exons_list:
            exon_name = exon_info[0]
            exon_last_info = exon_info[2]
            if exon_name not in exons_dic:
                exons_dic[exon_name] = []
            exons_dic[exon_name].append(
                (al_obj_in_hsp.position_query, al_obj_in_hsp.position_query_local, exon_last_info))

    # gff     
    for exon, gff_exon_info_list in exons_dic.iteritems():
        gff3_atributes = "GeneName=" + al_obj_in_hsp.subject_name + ";ExonName=" + exon + ";" + exon_last_info

        # global
        set_position_query = set(x[0] for x in gff_exon_info_list)
        save_exoninfo_in_gff(out_gff3_global_hlr,
                             al_obj_in_hsp.query_name,
                             gff3_source,
                             gff3_feature,
                             str(min(set_position_query)),
                             str(max(set_position_query)),
                             str(al_obj_in_hsp.score),
                             al_obj_in_hsp.strand,
                             gff3_frame,
                             gff3_atributes)

        # local
        set_position_query_local = set(x[1] for x in gff_exon_info_list)
        save_exoninfo_in_gff(out_gff3_local_hlr,
                             al_obj_in_hsp.query_geneName_local,
                             gff3_source,
                             gff3_feature,
                             str(min(set_position_query_local)),
                             str(max(set_position_query_local)),
                             str(al_obj_in_hsp.score),
                             al_obj_in_hsp.strand,
                             gff3_frame,
                             gff3_atributes)

    """ Output data model

        001    .    blast_hit_model    444    444    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000408762;Last=False
        001    .    blast_hit_model    445    491    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000385819;Last=False
        001    .    blast_hit_model    492    607    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000377032;Last=False
        001    .    blast_hit_model    608    714    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000373100;Last=False
        001    .    blast_hit_model    715    834    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000423072;Last=False
        001    .    blast_hit_model    835    930    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000377711;Last=False


        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    1    1    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000408762;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    2    48    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000385819;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    49    164    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000377032;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    165    271    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000373100;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    272    391    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000423072;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    392    487    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000377711;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    488    615    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000409016;Last=False
        001__q[651:4048]_s[444:3826]_ENSXETG00000026876    .    blast_hit_model    616    709    0.0    +    .    GeneName=ENSXETG00000026876;ExonName=ENSXETE00000414542;Last=False


        gff3 model:

        ENSXETG00000008118    .    gene_model_exons    1    33    .    -1    .    Name=ENSXETE00000362374;Last=False
        ENSXETG00000008118    .    gene_model_exons    34    135    .    -1    .    Name=ENSXETE00000100066;Last=False

        seqname - The name of the sequence. Typically a chromosome or a contig. Argo does not care what you put here. It will superimpose gff features on any sequence you like.
        source - The program that generated this feature. Argo displays the value of this field in the inspector but does not do anything special with it.
        feature - The name of this type of feature. The official GFF3 spec states that this should be a term from the SOFA ontology, but Argo does not do anything with this value except display it.
        start - The starting position of the feature in the sequence. The first base is numbered 1.
        end - The ending position of the feature (inclusive).
        score - A score between 0 and 1000. If there is no score value, enter ".".
        strand - Valid entries include '+', '-', or '.' (for don't know/don't care).
        frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'. Argo does not do anything with this field except display its value.
        GFF3: grouping attributes Attribute keys and values are separated by '=' signs. Values must be URI encoded.quoted. Attribute pairs are separated by semicolons. Certain, special attributes are used for grouping and identification (See below). This field is the one important difference between GFF flavors.

     """


def get_hsp_alignment_object_list(hsp, alignment_geneName, query_geneName, query_geneName_local):
    """ Generates list of hsp alignments. """

    position_subject = hsp.sbjct_start
    position_query = hsp.query_start  # starting position for query, WITHOUT gaps
    position_query_local = 1

    # test for strand
    gff3_strand = "+"
    if hsp.sbjct_start > hsp.sbjct_end:
        gff3_strand = "-"

    # storage for output
    alignments_objects_list = []

    for index in xrange(len(str(hsp.sbjct))):  # length of alignment, with gaps!

        letter_subject = str(hsp.sbjct)[index]
        letter_query = str(hsp.query)[index]

        # case 1:
        if letter_subject == "-" and letter_query == "-":
            print("WARNING: gap in sbj & query, it should not be like this!")

            # position_subject += 1
            continue

        # case 2
        elif letter_subject != "-" and letter_query != "-":  # this idx has exon:

            al_base = AlignmentBase(letter_subject, letter_query, position_subject, position_query,
                                    position_query_local, query_geneName, alignment_geneName, hsp.expect, gff3_strand,
                                    query_geneName_local)

            alignments_objects_list.append(al_base)

            # minus strand, thus base position countdown
            if gff3_strand == "-":
                position_subject -= 1
            else:
                position_subject += 1

            position_query += 1
            position_query_local += 1

            continue

        # case 3
        elif letter_subject == "-" and letter_query != "-":

            al_base = AlignmentBase(letter_subject, letter_query, position_subject, position_query,
                                    position_query_local, query_geneName, alignment_geneName, hsp.expect, gff3_strand,
                                    query_geneName_local)
            alignments_objects_list.append(al_base)

            position_query += 1
            position_query_local += 1

            continue

        # case 4
        elif letter_subject != "-" and letter_query == "-":

            # minus strand, thus base position countdown
            if gff3_strand == "-":
                position_subject -= 1
            else:
                position_subject += 1

            continue

    return alignments_objects_list


def main(args=[]):
    usage = '''
    usage: %prog [options] arg \nProgram parses blast XML file, translates exons boundaries (annotated on reference sequences) to query sequences (e.g. transcript)"
            '''

    parser = OptionParser(usage, version='%prog version 1.0')

    parser.add_option("-r", "--reference_fasta", dest="REFERENCE_FASTA", help="reference in fasta format")
    parser.add_option("-g", "--gff_reference_fasta", dest="GFF_REFERENCE_FASTA", help="annotation for reference in gff3 format")
    parser.add_option("-q", "--query_fasta", dest="QUERY_FASTA", help="query in fasta format")

    parser.add_option("-b", "--blast_db_path_and_name", dest="BLAST_DB_PATH_AND_NAME", help="blast+ database" ''', default="blast_out.xml"''')
    parser.add_option("-v", "--blast_xml_file", dest="BLAST_XML_FILE", help="blast results in xml file format", default="blast_out.xml")
    parser.add_option("-f", "--output_file", dest="OUTPUT_FILE", help="output file", action="store", type="string", default=str(__name__) + ".txt")
    parser.add_option("-o", "--output_folder", dest="OUTPUT_FOLDER", help="output folder", default="./")
    parser.add_option("-e", "--e_value_thresh", dest="E_VALUE_THRESH", help="threshold e-value", default=1e-8)
    parser.add_option("-a", "--only_best_Alignment", dest="ONLY_BEST_ALIGNMENT", help="take only 1, best q-s pair", default=True)

    parser.add_option("-w", "--blast_word_size", dest="BLAST_WORD_SIZE", help="blast word_size", default=11)
    parser.add_option("-n", "--blast_num_threads", dest="BLAST_NUM_THREADS", help="number of threads", default=2)

    parser.add_option("-m", "--blast_match_score", dest="BLAST_MATCH_SCORE", help="reward for nt match", default=1)
    parser.add_option("-s", "--blast_mismatch_score", dest="BLAST_MISMATCH_SCORE", help="penalty for nt mismatch", default=-3)
    parser.add_option("-y", "--blast_gap_open", dest="BLAST_GAP_OPEN", help="cost of opening a gap", default=5)
    parser.add_option("-x", "--blast_gap_extend", dest="BLAST_GAP_EXTEND", help="cost of gap extension", default=2)

    (options, arg) = parser.parse_args(args)

    # --- Entering program

    t_st = time.time()

    if not os.path.isdir(options.OUTPUT_FOLDER):
        sys.stdout.write('\nWrong output directory!')
        return

    os.chdir(options.OUTPUT_FOLDER)

    logging_file = "log_output"
    if options.OUTPUT_FILE != "":
        logging_file = options.OUTPUT_FILE
    log_info_hlr = open(options.OUTPUT_FOLDER + os.sep + logging_file + ".log", "w")

    log_info = "Entering program: {}\n".format(os.path.basename(__file__))
    sys.stdout.write(log_info)
    log_info_hlr.write(log_info)

    log_info = "\nUsed options: {}\n".format("\n".join(str(options).split(",")))
    sys.stdout.write(log_info)
    log_info_hlr.write(log_info)

    # --- workspace

    s = os.path.join(os.path.dirname(__file__), '.')
    os.chdir(s)
    print os.getcwd()

    # --- parsing gff3 file

    log_info = "Parsing {}\n".format(options.GFF_REFERENCE_FASTA)
    sys.stdout.write(log_info)
    log_info_hlr.write(log_info)

    gff_dic = {}
    s_prev_gen = ""
    gff_ref_hlr = open(options.GFF_REFERENCE_FASTA, "r")
    for line_gff in gff_ref_hlr:
        line_gff_list = line_gff.split("\t")

        gene_name = line_gff_list[0]
        exon_start = int(line_gff_list[3])
        exon_end = int(line_gff_list[4])
        exon_strand = line_gff_list[6]
        exon_name = line_gff_list[8].split(";")[0].split("=")[1]
        last_exon = line_gff_list[8].split(";")[1].strip()

        info_pack = [exon_name, exon_strand, last_exon]

        if gene_name not in gff_dic:
            coord_exons_dic = {}
            for x_coord in range(exon_start, exon_end + 1):
                coord_exons_dic[x_coord] = [info_pack]
            gff_dic[gene_name] = coord_exons_dic
        else:
            coord_exons_dic = gff_dic[gene_name]
            for x_coord in range(exon_start, exon_end + 1):
                if x_coord not in coord_exons_dic:
                    coord_exons_dic[x_coord] = []
                coord_exons_dic[x_coord].append(info_pack)

        if gene_name != s_prev_gen:
            s_prev_gen = gene_name
            log_info = "Parsing gff for gene: {}\n".format(s_prev_gen)
            sys.stdout.write(log_info)
            log_info_hlr.write(log_info)

    # --- blast analysis for each sequence

    from Bio.Blast.Applications import NcbiblastnCommandline

    blast_db_source = options.BLAST_DB_PATH_AND_NAME
    blastx_cline = NcbiblastnCommandline(query=options.QUERY_FASTA, db=blast_db_source,
                                         evalue=float(options.E_VALUE_THRESH), outfmt=5, out=options.BLAST_XML_FILE,
                                         word_size=options.BLAST_WORD_SIZE,
                                         num_threads=options.BLAST_NUM_THREADS, reward=options.BLAST_MATCH_SCORE,
                                         penalty=options.BLAST_MISMATCH_SCORE, gapopen=options.BLAST_GAP_OPEN,
                                         gapextend=options.BLAST_GAP_EXTEND)
    stdout, stderr = blastx_cline()

    # --- analysis of blast alignment

    out_file_core = "exons_alignment_by_blast_out"
    out_local_gff3_hlr = open(options.OUTPUT_FOLDER + os.sep + out_file_core + "_local.gff3", "w")
    out_global_gff3_hlr = open(options.OUTPUT_FOLDER + os.sep + out_file_core + "_global.gff3", "w")
    out_fasta_hlr = open(options.OUTPUT_FOLDER + os.sep + out_file_core + ".fasta", "w")

    result_handle = open(options.BLAST_XML_FILE)
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < float(options.E_VALUE_THRESH):

                    alignment_geneName = str(alignment.hit_def)

                    print('sequence:', alignment.title)

                    print("len hsp.query", len(hsp.query))
                    print("len hsp.sbjct", len(hsp.sbjct))
                    print("len hsp.match", len(hsp.match))

                    print("hsp.sbjct", str(hsp.sbjct))
                    print("hsp.match", str(hsp.match))
                    print("hsp.query", str(hsp.query))

                    # coordinates: subject
                    print("hsp.sbjct_start", hsp.sbjct_start)
                    print("hsp.sbjct_start", hsp.sbjct_end)

                    # coordinates: query
                    print("hsp.sbjct_start", hsp.query_start)
                    print("hsp.sbjct_start", hsp.query_end)

                    coord_exons_dic = gff_dic[alignment_geneName]

                    # generate alignment objects list
                    query_geneName = str(blast_record.query)
                    query_geneName_local = query_geneName + "__q[" + str(hsp.query_start) + ":" + str(hsp.query_end) + "]" + "_s[" + str(hsp.sbjct_start) + ":" + str(
                        hsp.sbjct_end) + "]" + "_" + alignment_geneName
                    alignment_object_list = get_hsp_alignment_object_list(hsp, alignment_geneName, query_geneName, query_geneName_local)

                    query_seq = "".join([xx.q for xx in alignment_object_list])
                    out_fasta_hlr.write(">" + query_geneName_local + "\n")
                    out_fasta_hlr.write(query_seq + "\n")

                    # set exons info into alignment objects
                    for al_obj_in_hsp in alignment_object_list:
                        if al_obj_in_hsp.position_subject in coord_exons_dic:
                            al_obj_in_hsp.set_exons(coord_exons_dic[al_obj_in_hsp.position_subject])

                    # global & local gff output
                    extract_and_write_gff(alignment_object_list, out_global_gff3_hlr, out_local_gff3_hlr)

            if options.ONLY_BEST_ALIGNMENT:
                break

    out_local_gff3_hlr.close()
    out_global_gff3_hlr.close()
    out_fasta_hlr.close()

    # --- closing program

    t_end = time.time()

    sys.stdout.write("\n\nWork done...")
    sys.stdout.write("\nProcess time [s]: " + str(t_end - t_st))


if __name__ == "__main__":
    main(sys.argv[1:])
