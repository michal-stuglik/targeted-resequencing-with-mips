'''
Created on Jul 3, 2014

@author: ania.fijarczyk
'''

from collections import defaultdict
import random


def uniq(inlist):  # wybor unikatowych elementow z listy redundantnych elementow
    # order preserving
    uniques = []
    for item in inlist:
        if item not in uniques:
            uniques.append(item)
    return uniques


def duplicates(inlist):
    dup = []
    for ele in inlist:
        if inlist.count(ele) > 1:
            dup.append(ele)
    return uniq(dup)


def non_dupl(inlist):
    nondup = []
    for ele in inlist:
        if inlist.count(ele) == 1:
            nondup.append(ele)
    return nondup


def merge_lists(big_list):
    new_list = []
    for ele in big_list:
        new_list += ele
    return new_list


def readGff(filename):
    fh = open(filename, 'r')
    linie = fh.readlines()
    k = [ele.split() for ele in linie]
    # in gff start is 1 (first base is 1), stop base inclusive
    info = [['\t'.join(ele) + '\n'] + [ele[0]] + [ele[8].split(';')[0].split('=')[1]] + [int(ele[3])] + [int(ele[4])] for ele in k]

    # create dictionary:
    # {(43, 'ENSXETG00000017426'): ['ENSXETG00000017426\t.\tgene_model_exons\t1\t144\t.\t1\t.\tName=ENSXETE00000392567;Last=False\n', 'ENSXETG00000017426', 'ENSXETE00000392567', 1, 145],  ... }
    i = 0
    D = {}
    for linijka in info:
        D[i, linijka[1]] = linijka
        i += 1
    return D


def filter_exons(big_list):
    # group exons into tuples according to repetitive start and stop coordinates and pick the longest one from the group
    rep_start_coord = duplicates([ele[1] for ele in big_list])
    rep_stop_coord = duplicates([ele[2] for ele in big_list])

    N = []
    for exon in big_list:
        if (exon[1] not in rep_start_coord) and (exon[2] not in rep_stop_coord):
            N.append(exon)

    T_start = []
    T_stop = []
    if rep_start_coord:
        for duplicate in rep_start_coord:
            t = [ele for ele in big_list if ele[1] == duplicate]
            maks_length = max([ele[3] for ele in t])
            the_longest_exons = [ele for ele in t if ele[3] == maks_length]
            one_longest_exon = random.sample(the_longest_exons, 1)[0]
            T_start.append(one_longest_exon)

    if rep_stop_coord:
        for duplicate in rep_stop_coord:
            t = [ele for ele in big_list if ele[2] == duplicate]
            maks_length = max([ele[3] for ele in t])
            the_longest = [ele for ele in t if ele[3] == maks_length]
            one_longest = random.sample(the_longest, 1)[0]
            T_stop.append(one_longest)



    # remove identical exons
    T_uniq = uniq(T_start + T_stop + N)
    if T_uniq:
        return T_uniq
    else:
        return big_list


def writeDic(small_dict):
    wh = open('filter_exons.out', 'w')
    for ele in sorted(small_dict.keys()):
        wh.write(small_dict[ele])
    wh.flush()
    wh.close()

    '''
    P = []
    T = []
    for gene in small_dict.keys():
        for ele in small_dict[gene]:
            p = gene, ele[0], ele[1], ele[2]
            P.append(p)
    My_table=sorted(P, key=lambda a: (a[0],a[2]))
    for ele in My_table:
        t = ele[0]+'\t'+ele[1]+'\t'+str(ele[2])+'\t'+str(ele[3])+'\n'
        T.append(t)
    T_print = ''.join(T)
    wh.write(T_print)
    '''


def filter_dictionary(main_dict):
    # create dictionary with gene ID as a key : exonID, start, stop, length:
    # defaultdict(<type 'list'>, {'ENSXETG00000014799': [['ENSE00002619574', 640, 737, 97], ['ENSE00002831632', 473, 576, 103], ['ENSE00002609838', 803, 887, 84], ...]}
    d = defaultdict(list)
    for ele in main_dict.keys():
        d[main_dict[ele][1]].append(main_dict[ele][2:] + [main_dict[ele][4] - main_dict[ele][3] + 1])

    R = {}
    for gene in d.keys():  # gene= [['ENSE00002619574', 640, 737, 97], ['ENSE00002831632', 473, 576, 103], ['ENSE00002609838', 803, 887, 84],...]

        # ref = d['ENSG00000005844']
        ref = d[gene]
        # check if there are duplicated start or stop coordinates
        starts = [ele[1] for ele in ref]
        stops = [ele[2] for ele in ref]
        dup_starts = duplicates(starts)
        dup_stops = duplicates(stops)


        # remove redundant exons (with identical left and right coordinates), ie pick one randomly
        f = defaultdict(list)
        for a, b, c, e in ref:
            f[(b, c, e)].append(a)
        f_sampled = {}
        for coordinates in f.keys():
            f_sampled[coordinates] = random.sample(f[coordinates], 1)
        F = []
        for ele in f_sampled.keys():
            F.append((f_sampled[ele] + [ele[0], ele[1], ele[2]]))

        # split dataset into those exons which have uniq coordinates and those with repetitive coordinates
        D = []
        S = []
        for exon in F:
            if (exon[1] in dup_starts) or (exon[2] in dup_stops):
                D.append(exon)
            elif (exon[1] not in dup_starts) and (exon[2] not in dup_stops):
                S.append(exon)


        # check if there are repetitive coordinates
        rep_starts = duplicates([ele[1] for ele in D])
        rep_stops = duplicates([ele[2] for ele in D])
        all_reps = rep_starts + rep_stops


        # if there are no exons with repetitive coordinates
        if all_reps == []:
            final_exons = D + S
        else:
            # loop for filtering exons
            while all_reps:
                D = filter_exons(D)

                new_rep_starts = duplicates([ele[1] for ele in D])
                new_rep_stops = duplicates([ele[2] for ele in D])
                all_reps = new_rep_starts + new_rep_stops

        final_exons = D + S

        R[gene] = final_exons

    # x = F['ENSG00000005844']
    return R


if __name__ == '__main__':
    import sys
    # fname = 'xenopus.gff3'
    # fname = 'joined_gene_model.gff3'
    fname = sys.argv[1]

    exon_dict = readGff(fname)
    L = {}
    Nums = {}
    for ele in exon_dict.keys():
        L[ele[0]] = exon_dict[ele][0]  # dictionary {1 : info}
        Nums[exon_dict[ele][2]] = ele[0]  # dictionary {exon : 1}

    f = filter_dictionary(exon_dict)
    filtered_exon_list = [ele[0] for ele in merge_lists(f.values())]

    # select f exons from exon_dict
    New_dict = {}
    for exon in filtered_exon_list:
        New_dict[Nums[exon]] = L[Nums[exon]]

    w = writeDic(New_dict)

    # print New_dict
