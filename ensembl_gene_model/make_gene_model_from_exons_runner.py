import os
import make_gene_model_from_exons


def main():
    s = r"../sample"

    args = [
        "-f", s + os.sep + "xenopus_exons.fasta",
        "-e", s + os.sep + "xenopus_exons_info.txt",
        "-o", s
    ]

    make_gene_model_from_exons.main(args)


if __name__ == "__main__":
    main()
