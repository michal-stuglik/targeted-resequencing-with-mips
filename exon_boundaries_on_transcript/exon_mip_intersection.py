'''
Created on Mar 2, 2015

@author: ania.fijarczyk


This script reads 'exons_alignment_by_blast_out_global.gff3' file, looks for overlapping exons in each reference,
and for each overlapping pair of exons, splits them into non-overlapping regions.

gff3 file (input):
001    .    megablast_hit    525    699    5e-88    +    .    GeneName=001;ExonName=Contig181_1
001    .    megablast_hit    6508    7527    0.0    +    .    GeneName=001;ExonName=Contig431_2
001    .    megablast_hit    4145    5588    0.0    +    .    GeneName=001;ExonName=Contig432_3

Reference is given in the first column, exon names are given in the last column ('ExonName')

USAGE: python exon_mip_intersection.py

The output is called 'exon_mip_intersection.out' and consists of three columns:
reference name, exon(region) start position, exon(region) stop position (numeration as in a gff file)
'''

from collections import defaultdict
import operator


def sort_table(table, col=0):
    return sorted(table, key=operator.itemgetter(col))


def addLists(lista):
    N = []
    for ele in lista:
        E = list(ele)
        N += E
    return N


def uniq(inlist):  # choose unique elements from a list
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques

def find_exon_name(instring):
    info_list = instring.split(';')
    info_two = [ele.split('=') for ele in info_list]
    D = {a:b for a,b in info_two}
    exon_name = D['ExonName']
    return exon_name


def readMe(filename):
    fh = open(filename, 'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    n = [(ele[0], int(ele[3]), int(ele[4]), find_exon_name(ele[8])) for ele in k]
    d = defaultdict(list)
    for hit in n:
        d[hit[0]].append(hit[1:])
    return d


def split_coordinates(lista):
    n = len(lista) - 1
    P = []
    for i in range(n):
        p = (lista[i], lista[i + 1] - 1)
        P.append(p)
    return P


if __name__ == '__main__':

    fname = 'exons_alignment_by_blast_out_global.gff3'
    # fname = sys.argv[1]
    r = readMe(fname)

    # list of merged contigs
    R = {}
    C = {}
    for ref in r.keys():
        Contigs = {nazwa: (start, stop) for start, stop, nazwa in r[ref]}
        C[ref] = Contigs
        old_unsorted = [(start, stop) for start, stop, nazwa in r[ref]]
        Old_list_whole = sort_table(old_unsorted, 0)

        Old_list = Old_list_whole[1:]
        New_list = [Old_list_whole[0]]

        num = len(Old_list)

        for exon in range(num):
            first = New_list[-1]
            second = Old_list[0]
            first_set = range(first[0], first[1] + 1)
            second_set = range(second[0], second[1] + 1)
            common = list(set(first_set) & set(second_set))
            common_sort = common.sort()
            if common:
                ab = uniq(first_set + second_set)
                ab_contig = min(ab), max(ab)

                del New_list[-1]
                New_list.append(ab_contig)

            else:
                New_list.append(second)

            Old_list.remove(second)

        R[ref] = New_list

    # for each merged contig put start and stop positions and split into adjecent intervals
    D = {contig: [range(ele[0], ele[1] + 1) for ele in R[contig]] for contig in R.keys()}

    G = {}
    for gen in D.keys():
        koordynaty = [(k[0], k[1]) for k in r[gen]]
        S = []
        for merged_contig in D[gen]:
            s = [(coord[0], coord[1] + 1) for coord in koordynaty if (coord[0] in merged_contig) and (coord[1] in merged_contig)]
            sum_sort_coord = uniq(sorted(addLists(s)))
            split_regions = split_coordinates(sum_sort_coord)
            S.append(split_regions)

        G[gen] = addLists(S)

    G_formated = {g: [(str(ele[0]) + '\t' + str(ele[1])) for ele in G[g]] for g in G.keys()}

    wh = open('exon_mip_intersection.out', 'w')
    for reference in sorted(G_formated.keys()):
        for inter in G_formated[reference]:
            wh.write(reference + '\t' + inter + '\n')

    wh.flush()
    wh.close()

    # print G_formated

