'''
Created on Jun 30, 2014

@author: ania.fijarczyk
'''
from collections import defaultdict
import itertools
from sets import Set


def merge_lists(big_list):
    new_list = []
    for ele in big_list:
        new_list += ele
    return new_list


def sort_columns(lista):
    My_table = sorted(lista, key=lambda a: (a[2]))
    return My_table


def readMe(filename):
    fh = open(filename, 'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    # info for .gff3, where exon ID in the second position
    info = [(ele[0], ele[3], ele[4], ele[8].split(';')[1].split('=')[1]) for ele in k]
    d = defaultdict(list)
    for ele in info:
        d[ele[0]].append(ele[1:])

    return d


def write_dict_sorted(small_dict, ref):
    P = []
    for ele in small_dict.keys():
        p = ref, ele, small_dict[ele][0], small_dict[ele][1]
        P.append(p)
    P_sort = sort_columns(P)
    P_sort_format = [(ele[0] + '\t' + ele[1] + '\t' + str(ele[2]) + '\t' + str(ele[3]) + '\n') for ele in P_sort]
    return P_sort_format


def writeBed(big_dict):
    B = []
    for ele in sorted(big_dict.keys()):
        p = write_dict_sorted(big_dict[ele], ele)
        B += p
    B_to_write = ''.join(B)
    wh = open('exon.bed', 'w')
    # wh.write('reference\tEnsemble_exon_ID\texon_start\texon_stop\n')
    wh.write(B_to_write)
    wh.flush()
    wh.close()


def writeStat(small_dict):
    wh = open('exon_stat.out', 'w')
    for ele in sorted(small_dict.keys()):
        wh.write(ele + '\t' + str(small_dict[ele][0]) + '\t' + str(small_dict[ele][1]) + '\n')
    wh.flush()
    wh.close()


def writeOverlap(small_dict):
    wh = open('exon_overlap.out', 'w')
    for ele in sorted(small_dict.keys()):
        if small_dict[ele]:
            for pair in small_dict[ele]:
                wh.write(ele + '\t' + pair[0] + '\t' + pair[1] + '\t' + str(pair[2]) + '\t' + str(pair[3]) + '\t' + str(pair[4]) + '\t' + str(pair[5]) + '\t' + str(pair[6]) + '\n')
    wh.flush()
    wh.close()


class EXON:
    def __init__(self, lista):
        self.info = lista
        self.exons = {}
        self.bed = {}
        self.lengths = {}
        self.common = []
        self.pairs = []

        self.get_exons(lista)
        self.get_common(self.exons)

    def get_exons(self, start_stop_exon_lista):
        D = defaultdict(list)
        for ele in start_stop_exon_lista:
            D[ele[2]].append((int(ele[0]), int(ele[1])))
        Range = {}
        Bed = {}
        Lengths = {}
        for ele in D.keys():
            t = merge_lists(D[ele])
            Range[ele] = range(min(t), max(t) + 1)
            Bed[ele] = min(t) - 1, max(t) - 1
            Lengths[ele] = len(range(min(t), max(t) + 1)) - 1
        self.exons = Range
        self.bed = Bed
        self.lengths = Lengths

    def get_common(self, leks):
        exon_ids = leks.keys()
        pary = itertools.combinations(exon_ids, 2)
        p = list(pary)
        R = []
        for dwojka in p:
            first_set = Set(leks[dwojka[0]])
            second_set = Set(leks[dwojka[1]])
            common_part = first_set & second_set
            if len(common_part):
                # r = dwojka[0], dwojka[1], len(common_part)
                r = dwojka[0], dwojka[1], min(leks[dwojka[0]]), max(leks[dwojka[0]]), min(leks[dwojka[1]]), max(leks[dwojka[1]]), len(common_part)
                R.append(r)
        self.pairs = R


if __name__ == '__main__':
    import sys

    fname = sys.argv[1]
    # fname = 'exons_alignment_by_blast_out_global.gff3'


    r = readMe(fname)
    R = {}
    Bed = {}
    L = {}
    for ref in r.keys():
        exon = EXON(r[ref])
        pary = exon.pairs
        R[ref] = pary  # overlapping exons
        Bed[ref] = exon.bed
        dlugosci = exon.lengths.values()
        L[ref] = len(dlugosci), sum(dlugosci)
    # w1 = writeBed(Bed) # coordinates have -1 positions
    w2 = writeStat(L)
    w3 = writeOverlap(R)
    # print R
