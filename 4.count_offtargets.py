from sys import argv
from collections import defaultdict

FILE_CNT = 21

class Info_Transcript:
    def __init__(self, pk, cds_pos):
        self.pk = pk
        self.cds_pos = [[int(i) for i in posstr.split(',')] for posstr in cds_pos]

def find_transcript(transcripts, pos, direction, seqlen=23):
    if direction == '+':
        bp = pos + seqlen - 6
    else:
        bp = pos + 6
    for transc in transcripts:
        for cds_start, cds_end in transc.cds_pos:
            if cds_start < bp < cds_end:
                return transc.pk
    return -1

def main():
    cnt_dic = defaultdict(lambda: [0, 0, 0])
    p = argv[1].rstrip('/')

    print('Reading required files...')
    seq_pk_dic = {}
    with open(p+'/info_sequences.txt') as f:
        for line in f:
            pk, seq = line.strip().split('\t')[:2]
            seq_pk_dic[seq] = pk

    transc_per_chrom_dic = defaultdict(lambda: [])
    with open(p+'/info_transcripts.txt') as f:
        for line in f:
            entries = line.strip().split('\t')
            info_transc = Info_Transcript(entries[0], entries[4:])
            transc_per_chrom_dic[entries[2]].append(info_transc)

    pks = [1,1,1]
    info_files = [ open(p+'/info_offtargets_0.txt', 'w'),
                   open(p+'/info_offtargets_1.txt', 'w'),
                   open(p+'/info_offtargets_2.txt', 'w') ]

    print('Counting off-targets and locating its origin...')
    for i in range(FILE_CNT):
        with open(p+'/outs_%d.txt'%(i+1)) as f, open(p+'/offtarget_notfound.txt', 'w') as fo:
            for line in f:
                entries = line.strip().split('\t')
                num_off = int(entries[5])
                try:
                    info_files[num_off].write('\t'.join([str(pks[num_off]), str(seq_pk_dic[entries[0]])] + entries[1:-1]) + '\t' + str(find_transcript(transc_per_chrom_dic[entries[1]], int(entries[2]), entries[4])) + '\n')
                    cnt_dic[entries[0]][num_off] += 1
                    pks[num_off] += 1
                except KeyError:
                    fo.write(entries[0] + '\n')

    for f in info_files:
        f.close()

    print('Printing results...')
    with open(p+'/info_offtargets.txt', 'w') as fo:
        for target in cnt_dic:
            fo.write(str(seq_pk_dic[target]) + '\t' + '\t'.join([str(i) for i in cnt_dic[target]]) + '\n')

if __name__ == '__main__':
    main()
