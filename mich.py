from re import findall
from math import exp

try:
    from string import maketrans
    t = maketrans("ATGC", "TACG")
except:
    t = str.maketrans("ATGC", "TACG")

import re

def finditer_everything(pattern, string, flags=0):
    pattern = re.compile(pattern, flags=flags)
    pos = 0
    m = pattern.search(string, pos)
    while m is not None:
        yield m
        pos = m.start() + 1
        m = pattern.search(string, pos)


def target_yield(exon_seq, n=60, seed_len=20, pam_seq="NGG"):
    if n%2 == 1:
        raise Exception
    exon_seq = str(exon_seq)
    pam_seq = str(pam_seq)
    pattern = seed_len*'N' + pam_seq
    pattern = "(" + pattern + ")|(" + pattern.translate(t)[::-1] + ")"
    pattern = pattern.replace('N', '[AGTC]').replace('R', '[AG]').replace('W', '[AT]').replace('M', '[AC]')

    for m in finditer_everything(pattern, exon_seq):
        for i in (1, 2,):
            if m.group(i):
                if i == 1:
                    seq_start = m.start()+seed_len-3-n/2
                    seq_RGEN = m.group(i)
                    direction = '+'
                else:
                    seq_start = m.start()+len(pam_seq)+3-n/2
                    seq_RGEN = m.group(i).translate(t)[::-1]
                    direction = '-'

                seq_end = seq_start + n
                seq_long_pre, seq_long_post = '', ''

                if seq_start < 0:
                    seq_long_pre = '-'*abs(seq_start)
                    seq_start = 0
                if seq_end > len(exon_seq):
                    seq_long_post = '-'*(seq_end-len(exon_seq))
                    seq_end = len(exon_seq)

                seq_long = seq_long_pre + exon_seq[seq_start:seq_end] + seq_long_post
                yield (m.group(i), direction, seq_RGEN, m.start(), seq_long)

def mich_yield(mich_seq, left, length_weight):
    #if mich_seq[0] == '-' or mich_seq[-1] == '-':
    #    yield 0, 0
    sum_score_3=0
    sum_score_not_3=0

    right=len(mich_seq)-int(left)

    dup_list = []
    for k in range(2,left)[::-1]:
        for j in range(left,left+right-k+1):
            for i in range(0,left-k+1):
                if mich_seq[i:i+k]==mich_seq[j:j+k]:
                    length=j-i
                    dup_list.append( (mich_seq[i:i+k], i, i+k, j, j+k, length) )

    for i, dup in enumerate(dup_list):
        n=0
        score_3=0
        score_not_3=0
        
        scrap=dup[0]
        left_start=dup[1]
        left_end=dup[2]
        right_start=dup[3]
        right_end=dup[4]
        length=dup[5]

        for j in range(i):
            left_start_ref=dup_list[j][1]
            left_end_ref=dup_list[j][2]
            right_start_ref=dup_list[j][3]
            right_end_ref=dup_list[j][4]

            if (left_start >= left_start_ref) and (left_end <= left_end_ref) and (right_start >= right_start_ref) and (right_end <= right_end_ref):
                if (left_start - left_start_ref)==(right_start - right_start_ref) and (left_end - left_end_ref)==(right_end - right_end_ref):
                    n+=1
            else: pass

        if n == 0:
            length_factor = round(1/exp((length)/(length_weight)),3)
            num_GC=len(findall('G',scrap))+len(findall('C',scrap))
            score=100*length_factor*((len(scrap)-num_GC)+(num_GC*2))
            if (length % 3)==0:
                flag_3 = 0
            elif (length % 3)!=0:
                flag_3 = 1

            yield ((mich_seq[0:left_end]+'-'*length+mich_seq[right_end:], scrap, str(length), 100*length_factor*((len(scrap)-num_GC)+(num_GC*2))), flag_3)

def calc_mich_score(mich_seq):
    mich_seq = mich_seq.upper().strip()

    if mich_seq[0] == '-' or mich_seq[-1] == '-':
        return [], "", ""
    tot_score = 0
    tot_not_score = 0
    tot_list = []
    for tup, flag_3 in mich_yield(mich_seq, int(len(mich_seq)/2), 20.0):
        tot_list.append(tup)
        if flag_3 == 0:
            tot_score += tup[3]
        else:
            tot_not_score += tup[3]
    mich_score = tot_score+tot_not_score
    if mich_score != 0:
        oof_score = str(tot_not_score*100.0/mich_score)
        tot_list.sort(key=lambda e: e[3], reverse=True)
        tot_list.insert(0, (mich_seq, "", 0, 0))
    else:
        oof_score = "NaN"
    return tot_list, str(mich_score), oof_score
