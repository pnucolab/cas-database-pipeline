import os
from sys import argv
import json

p = argv[1].rstrip('/')
tablename = argv[2]

def load_file(filename, info_dic):
    with open(p+'/'+filename) as f:
        for cnt, line in enumerate(f, start=1):
            entries = line.strip().split('\t')
            info_dic[int(entries[0])] = entries[1:]

class fixture:
    def __init__(self, p):
        if not os.path.isdir(p+"/fixtures"):
            os.mkdir(p+"/fixtures")
        self.fixture_item_count = 0
        self.fixture_count = 1
        self.fixture_file = open(p+"/fixtures/fixture_1.json", 'w')
        self.fixture_file.write('[')

    def __enter__(self):
        return self

    def append(self, s):
        if self.fixture_item_count == 1000000:
            self.fixture_file.write(']')
            self.fixture_file.close()
            self.fixture_count += 1
            self.fixture_file = open(p+"/fixtures/fixture_%d.json"%self.fixture_count, 'w')
            self.fixture_file.write('[')
            self.fixture_item_count = 0
        else:
            if self.fixture_item_count != 0:
                self.fixture_file.write(', ')
        self.fixture_file.write(s)
        self.fixture_item_count += 1

    def close(self):
        self.fixture_file.close()

    def __exit__(self, t, v, tr):
        try:
            self.fixture_file.close()
        except:
            pass
"""
def fixture_append(s):
    if fixture_item_count == 100000 or fixture_file == None:
        if fixture_file:
            fixture_file.write(']')
            fixture_file.close()
        fixture_file = open(p+"/fixture_%d.json"%fixture_count)
        fixture_file.write('[')
        fixture_count += 1
        fixture_item_count = 0
    else:
        fixture_file.write(', ')
    fixture_file.write(s)
    fixture_item_count += 1
"""


print("Loading all required data files...")
#info_genes = {}
#info_transcripts = {}
#info_targets = {}
#info_target_transcripts = {}
#info_sequences = {}
info_offtargets = {}
#info_offtargets_0 = {}
#info_offtargets_1 = {}
#info_offtargets_2 = {}

#load_file('info_genes.txt', info_genes)
#load_file('info_transcripts.txt', info_transcripts)
#load_file('info_targets.txt', info_targets)
#load_file('info_target_transcripts.txt', info_target_transcripts)
#load_file('info_sequences.txt', info_sequences)

load_file('info_cpp_offtargets.txt', info_offtargets)
#load_file('info_cpp_offtargets_0.txt', info_offtargets_0)
#load_file('info_cpp_offtargets_1.txt', info_offtargets_1)
#load_file('info_cpp_offtargets_2.txt', info_offtargets_2)

with fixture(p) as f:
    print("Generating fixtures (Genes)...")
    with open(p+"/info_genes.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, gene = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_gene"%tablename,
                                    "pk": pk,
                                    "fields": { "geneid": gene[0],
                                                "sym": gene[1],
                                                "desc": gene[2] } } ) )

    print("Generating fixtures (Transcripts)...")
    with open(p+"/info_transcripts.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, transc = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_transcript"%tablename,
                                    "pk": pk,
                                    "fields": { "transcriptid": transc[0],
                                                "chrom": transc[1],
                                                "gene": int(transc[2]) } } ) )

    print("Generating fixtures (Targets)...")
    with open(p+"/info_targets.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, tgt = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_target"%tablename,
                                    "pk": pk,
                                    "fields": { "pam": tgt[0],
                                                "chrom": tgt[1],
                                                "pos": int(tgt[2]),
                                                "dir": tgt[3],
                                                "oof_score": float(tgt[4]),
                                                "cov": int(tgt[5]),
                                                "transcripts": [int(i) for i in tgt[6:]] } } ) )

    print("Generating fixtures (Target-Transcripts)...")
    with open(p+"/info_target_transcripts.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, tgt_transc = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_target_transcript"%tablename,
                                    "pk": pk,
                                    "fields": { "target": int(tgt_transc[0]),
                                                "transcript": int(tgt_transc[1]),
                                                "cds_pos": float(tgt_transc[2])} } ) )

    print("Generating fixtures (Sequences)...")
    with open(p+"/info_sequences.txt") as fi, open(p+"/targets_notfound.txt", "w") as fo:
        fo.write('human_hg38\n'+'N'*21+'RG\n')
        for line in fi:
            entries = line.strip().split('\t')
            pk, seq = int(entries[0]), entries[1:]
            try:
                off = info_offtargets[pk]
            except:
                fo.write(seq[0] + ' 2\n')
            f.append( json.dumps( { "model": "cas_database.%s_sequence"%tablename,
                                    "pk": pk,
                                    "fields": { "seq": seq[0],
                                                "gc": float(seq[1]),
                                                "off_0": int(off[0]),
                                                "off_1": int(off[1]),
                                                "off_2": int(off[2]),
                                                "targets": [int(i) for i in seq[2:]] } } ) )

    print("Generating fixtures (Offtargets with 0)...")
    with open(p+"/info_cpp_offtargets_0.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, off = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_off0"%tablename,
                                    "pk": pk,
                                    "fields": { "sequence": int(off[0]),
                                                "chrom": off[1],
                                                "pos": int(off[2]),
                                                "seq": off[3],
                                                "dir": off[4],
                                                "transcript": int(off[5]) } } ) )

    print("Generating fixtures (Offtargets with 1)...")
    with open(p+"/info_cpp_offtargets_1.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, off = int(entries[0]), entries[1:] 
            f.append( json.dumps( { "model": "cas_database.%s_off1"%tablename,
                                    "pk": pk,
                                    "fields": { "sequence": int(off[0]),
                                                "chrom": off[1],
                                                "pos": int(off[2]),
                                                "seq": off[3],
                                                "dir": off[4],
                                                "transcript": int(off[5]) } } ) )

    print("Generating fixtures (Offtargets with 2)...")
    with open(p+"/info_cpp_offtargets_2.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, off = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.%s_off2"%tablename,
                                    "pk": pk,
                                    "fields": { "sequence": int(off[0]),
                                                "chrom": off[1],
                                                "pos": int(off[2]),
                                                "seq": off[3],
                                                "dir": off[4],
                                                "transcript": int(off[5]) } } ) )
