# # 1. download_gene_info_from_ensembl.py
from tqdm.notebook import tqdm
import pandas as pd

df = pd.read_csv('/home/hdh1028/cas/second.csv')

with open('/home/hdh1028/cas/SLmic1.0_gene_models.gff') as fp:
    data = fp.readlines()
    data = [item.replace('\n','') for item in data]
    for i in tqdm(range(len(data))):
        data[i] = data[i].split()
cds_dict = {}
for i in range(len(data)):
    #print(data[i][0])
    cds_dict[f'{data[i][0]}'] = []
# print(cds_dict)
for i in tqdm(range(len(data))):
    if data[i][2] == 'CDS':
        cds_dict[f'{data[i][0]}'].append([data[i][-1][3:13],int(data[i][3])-1,int(data[i][4])-1])#gene_name, 
cds_dict = {key:value for key,value in cds_dict.items() if key == 'SLmic1.0_chr1'
            or key == 'SLmic1.0_chr2'
            or key == 'SLmic1.0_chr3'
            or key == 'SLmic1.0_chr4'
            or key == 'SLmic1.0_chr5'
            or key == 'SLmic1.0_chr6'
            or key == 'SLmic1.0_chr7'
            or key == 'SLmic1.0_chr8'
            or key == 'SLmic1.0_chr9'
            or key == 'SLmic1.0_chr10' 
            or key == 'SLmic1.0_chr11'
            or key == 'SLmic1.0_chr12'}
cds_dict.keys()

f = open("gene_data", 'w')
f.write('#ensembl_gene_id\tensembl_transcript_id\texternal_gene_id\tchromosome_name\tdescription\tgenomic_coding_start\tgenomic_coding_end\n')
for i in cds_dict:
    for j in cds_dict[i]:
        ensembl_gene_id = j[0]
        ensembl_transcript_id = ''
        external_gene_id = ''
        chromosome_name = f'{i}'
        description = ''
        genomic_coding_start = j[1]
        genomic_coding_end = j[2]
        
        f.write(f'{ensembl_gene_id}\t{ensembl_transcript_id}\t{external_gene_id}\t{chromosome_name}\t{description}\t{genomic_coding_start}\t{genomic_coding_end}\n')
f.close()

f = open("/home/hdh1028/cas/cas_scripts/chrom_table", 'w')
f.write('#name_in_gene_data\tname_in_fasta\n')
for i in range(12):
    f.write(f'SLmic1.0_chr{i+1}\tSLmic1.0_chr{i+1}\n')
f.close()


# # 2. generate_targets.py
from sys import argv
from os.path import isfile
import re
import glob
from Bio import SeqIO
from Bio.SeqUtils import GC
import math
from collections import defaultdict
from mich import calc_mich_score


TARGET_NUM = 21
print(argv)
wdir = '/home/hdh1028/cas/cas_scripts' # Working directory 1 -> 0으로 수정
cdir = '/home/hdh1028/cas/SLmic1.0' # A directory containing fasta file



class Info_Guide:
    def __init__(self):
        self.pk = -1
        self.seq = ''
        self.gc = 0

class Info_Transcript:
    def __init__(self):
        self.pk = -1
        self.gene = 0
        self.transcid = ''
        self.chrname = ''

        # Internally used
        self.cds_positions = []
        self.targets = []
        self.target_transcripts = []

class Info_CDSRange:
    def __init__(self):
        self.pk = -1
        self.transcript = 0
        self.start = 0
        self.end = 0

class Info_Gene:
    def __init__(self):
        self.pk = -1
        self.geneid = ''
        self.name = ''
        self.desc = ''

        # Internelly used
        self.transcripts = []

class Info_Target:
    def __init__(self):
        self.pk = -1
        self.guide = 0
        self.pam = ''
        self.chrname = ''
        self.pos = 0
        self.direction = ''
        self.oof_score = 0
        self.coverage = 0
        self.transcripts = []

        # Internally used
        self.seq = ''

class Info_Target_Transcript:
    def __init__(self):
        self.pk = -1
        self.target = 0
        self.transcript = 0
        self.cds_pos = 0

def finditer_everything(pattern, string, flags=0):
    pattern = re.compile(pattern, flags=flags)
    pos = 0
    m = pattern.search(string, pos)
    while m is not None:
        yield m
        pos = m.start() + 1
        m = pattern.search(string, pos)

def findall_everything(pattern, string, flags=0):
    return [i.start() for i in finditer_everything(pattern, string, flags)]

def get_chrom_name(ng, chrom_table=None, chrom_rules=None):
    if chrom_table:
        return chrom_table[ng]
    else:
        for chrom_rule in chrom_rules:
            p = re.compile(chrom_rule[0])
            m = p.search(ng)
            if m and len(m.groups()) != 0:
                groups = []
                for group in m.groups():
                    if group.isdigit():
                        groups.append(int(group))
                    else:
                        groups.append(group)
                return chrom_rule[1].format(*groups)
            else:
                if ng == chrom_rule[0]:
                    return chrom_rule[1]
        return ng


def main():
    if not isfile(wdir+'/chrom_table') and not isfile(wdir+'/chrom_rule'):
        print("Couldn't find chromosome table, '"+wdir+"/chrom_table'.\n\nPlease create chromosome table as below:\n#name_in_gene_data\tname_in_fasta\n1\tchr1\n2\tchr2\n...\n\nOr you can define '"+wdir+"/chrom_rule' as below:\n#regex_name_in_gene_data\tformatted_name_in_fasta\n(.*)\tchr{0}\n\nPlease note that the regex pattern should contain at least one subgroup.")
        return
    chrom_table = {}
    chrom_rules = []
    if isfile(wdir+'/chrom_table'):
        with open(wdir+'/chrom_table') as f:
            for line in f:
                if line[0] == '#':
                    continue
                ng, nf = line.strip().split('\t')[:2]
                chrom_table[ng] = nf
    else:
        with open(wdir+'/chrom_rule') as f:
            for line in f:
                if line[0] == '#':
                    continue
                chrom_rules.append(line.strip().split('\t')[:2])

    chrom_data = {}
    cur_chrom = ''
    for fn in glob.glob(cdir+'/*'):
        print('Loading %s...'%fn)
        for seq_record in SeqIO.parse(fn, "fasta"):
            chrom_data[seq_record.id] = seq_record.seq

    print('Finding RGEN targets (NGG and CCN)...')
    err_chroms = []
    targets = []

    info_genes = defaultdict(lambda: Info_Gene())
    info_transcripts = defaultdict(lambda: Info_Transcript())
    info_targets = defaultdict(lambda: Info_Target())
    info_guides = defaultdict(lambda: Info_Guide())
    info_sequences = defaultdict(lambda: Info_Sequence())
    info_target_transcripts = defaultdict(lambda: Info_Target_Transcript())

    pk_gene = 1
    pk_transc = 1
    pk_tgt = 1
    pk_guide = 1
    pk_tgt_transc = 1

    with open(wdir+'/gene_data') as f:
        for linecnt, line in enumerate(f):
            if linecnt % 1000 == 0:
                print("Processing line #%d..."%(linecnt+1))
            if line[0] == '#':
                continue
            geneid, transcid, genename, chrname, desc, cds_start, cds_end = line.strip().split('\t')
            if chrom_table:
                chrname = get_chrom_name(chrname, chrom_table=chrom_table)
            else:
                chrname = get_chrom_name(chrname, chrom_rules=chrom_rules)
            if chrname in err_chroms:
                continue
            cds_start = int(cds_start)-1
            cds_end = int(cds_end)
            try:
                cds_seq = str(chrom_data[chrname][cds_start:cds_end]).upper()
            except KeyError:
                if not chrname in err_chroms:
                    print("Couldn't find %s in fasta files! Skipping..."%chrname)
                    err_chroms.append(chrname)
                continue

            info_gene = info_genes[geneid]
            if info_gene.pk == -1:
                info_gene.pk = pk_gene
                info_gene.geneid = geneid
                info_gene.name = genename
                info_gene.desc = desc
                pk_gene += 1

            info_transc = info_transcripts[transcid]
            if info_transc.pk == -1:
                info_transc.pk = pk_transc
                info_transc.transcid = transcid
                info_transc.chrname = chrname
                info_transc.gene = info_gene.pk
                info_gene.transcripts.append(transcid) # Assuming no duplication
                pk_transc += 1

            if not (cds_start, cds_end) in info_transc.cds_positions:
                info_transc.cds_positions.append( (cds_start, cds_end) ) # Internally Used

            seqs = []
            for pos in findall_everything(r'[AGTC]{5}[G]{2}', cds_seq):
                i = cds_start + pos - 16
                found_seq = str(chrom_data[chrname][i:i+23]).upper()
                seqs.append( (found_seq, i, '+') )

            for pos in findall_everything(r'[C]{2}[AGTC]{5}', cds_seq):
                i = cds_start + pos
                found_seq = str(chrom_data[chrname][i:i+23].reverse_complement()).upper()
                seqs.append( (found_seq, i, '-') )

            for seq, tgt_pos, direction in seqs:
                if not 'N' in seq:
                    pam = seq[-3:]
                    seq = seq[:-3]
                    targets.append(seq + "NNN 2")

                    info_guide = info_guides[seq]
                    if info_guide.pk == -1:
                        info_guide.pk = pk_guide
                        info_guide.seq = seq
                        info_guide.gc = GC(seq)
                        pk_guide += 1

                    tgt_id = ','.join( [chrname,str(tgt_pos),direction] )
                    info_tgt = info_targets[tgt_id]
                    if info_tgt.pk == -1:
                        info_tgt.pk = pk_tgt
                        info_tgt.guide = info_guide.pk
                        info_tgt.pam = pam
                        info_tgt.chrname = chrname
                        info_tgt.pos = tgt_pos
                        info_tgt.direction = direction
                        if direction == "+":
                            bp = tgt_pos + len(seq) - 3 # right side of bp
                        else:
                            bp = tgt_pos + 5 # left side of bp

                        if bp-30 < 0 or bp+30 > len(chrom_data[chrname]):
                            info_tgt.oof_score = -1
                        else:
                            info_tgt.oof_score = calc_mich_score(str(chrom_data[chrname][bp-30:bp+30]))[2]
                        pk_tgt += 1

                    if not tgt_id in info_transc.targets:
                        info_transc.targets.append(tgt_id)

                    if not info_transc.pk in info_tgt.transcripts:
                        info_tgt.transcripts.append(info_transc.pk)

    print('Calculating CDS position and coverage of each target...')
    for i in info_genes:
        for j in info_genes[i].transcripts:
            info_transc = info_transcripts[j]
            tot_cds_length = 0
            for cds_start, cds_end in info_transc.cds_positions:
                tot_cds_length += cds_end - cds_start

            for k in info_transc.targets:
                info_tgt = info_targets[k]
                info_tgt_transc = info_target_transcripts[str(info_tgt.pk)+','+str(info_transc.pk)]
                info_tgt_transc.pk = pk_tgt_transc
                info_tgt_transc.target = info_tgt.pk
                info_tgt_transc.transcript = info_transc.pk
                pk_tgt_transc += 1
                cds_length = 0
                info_transc.cds_positions.sort(key=lambda e: e[0])
                for cds_start, cds_end in info_transc.cds_positions:
                    if info_tgt.direction == '+':
                        bp = info_tgt.pos + 20 - 3 # right side of bp
                    else:
                        bp = info_tgt.pos + 6 # right side of bp

                    if cds_start < bp < cds_end:
                        cds_length += bp - cds_start
                        info_tgt_transc.cds_pos = (cds_length-1) / float(tot_cds_length) * 100.0 # calculated position cannot be the last seq position, thus cds_pos cannot be 100%.
                        info_tgt.coverage += 1
                        break
                    else:
                        cds_length += cds_end - cds_start

    print('Removing duplicate targets...')
    targets = list(set(targets))

    print('Writing to file...')
    targets_per_file = int(math.ceil(len(targets) / float(TARGET_NUM)))
    for cnt, i in enumerate(range(0, len(targets), targets_per_file), 1):
        with open(wdir+'/targets_'+str(cnt)+'.txt', 'w') as fo:
            fo.write(cdir+'\n')
            fo.write('N'*21 + 'RG\n')
            fo.write('\n'.join(targets[i:i+targets_per_file]))

    with open(wdir+'/info_genes.txt', 'w') as fo:
        for pk in info_genes:
            gene = info_genes[pk]
            fo.write('\t'.join( [str(gene.pk), gene.geneid, gene.name, gene.desc] ) + "\n" )
    with open(wdir+'/info_transcripts.txt', 'w') as fo:
        for pk in info_transcripts:
            transc = info_transcripts[pk]
            fo.write('\t'.join( [str(transc.pk), transc.transcid, transc.chrname, str(transc.gene)] ) + '\t' + '\t'.join( [str(s)+','+str(e) for s,e in transc.cds_positions] ) + "\n" )
    with open(wdir+'/info_targets.txt', 'w') as fo:
        for pk in info_targets:
            tgt = info_targets[pk]
            fo.write('\t'.join( [str(tgt.pk), str(tgt.guide), tgt.pam, tgt.chrname, str(tgt.pos), tgt.direction, str(tgt.oof_score), str(tgt.coverage)] ) + '\t' + '\t'.join([str(i) for i in tgt.transcripts]) + '\n')
    with open(wdir+'/info_guides.txt', 'w') as fo:
        for pk in info_guides:
            guide = info_guides[pk]
            fo.write('\t'.join( [str(guide.pk), guide.seq, str(guide.gc)] ) + '\n')
    with open(wdir+'/info_target_transcripts.txt', 'w') as fo:
        for pk in info_target_transcripts:
            tgt_transc = info_target_transcripts[pk]
            fo.write('\t'.join( [str(tgt_transc.pk), str(tgt_transc.target), str(tgt_transc.transcript), str(tgt_transc.cds_pos)] ) + '\n')


if not isfile(wdir+'/chrom_table') and not isfile(wdir+'/chrom_rule'):
    print("Couldn't find chromosome table, '"+wdir+"/chrom_table'.\n\nPlease create chromosome table as below:\n#name_in_gene_data\tname_in_fasta\n1\tchr1\n2\tchr2\n...\n\nOr you can define '"+wdir+"/chrom_rule' as below:\n#regex_name_in_gene_data\tformatted_name_in_fasta\n(.*)\tchr{0}\n\nPlease note that the regex pattern should contain at least one subgroup.")
chrom_table = {}
chrom_rules = []
if isfile(wdir+'/chrom_table'):
    with open(wdir+'/chrom_table') as f:
        for line in f:
            if line[0] == '#':
                continue
            ng, nf = line.strip().split('\t')[:2]
            chrom_table[ng] = nf
else:
    with open(wdir+'/chrom_rule') as f:
        for line in f:
            if line[0] == '#':
                continue
            chrom_rules.append(line.strip().split('\t')[:2])



chrom_data = {}
cur_chrom = ''
for fn in tqdm(glob.glob(cdir+'/*')):
    print('Loading %s...'%fn)
    for seq_record in SeqIO.parse(fn, "fasta"):
        chrom_data[seq_record.id] = seq_record.seq

print('Finding RGEN targets (NGG and CCN)...')
err_chroms = []
targets = []

info_genes = defaultdict(lambda: Info_Gene())
info_transcripts = defaultdict(lambda: Info_Transcript())
info_targets = defaultdict(lambda: Info_Target())
info_guides = defaultdict(lambda: Info_Guide())
info_sequences = defaultdict(lambda: Info_Sequence())
info_target_transcripts = defaultdict(lambda: Info_Target_Transcript())

pk_gene = 1
pk_transc = 1
pk_tgt = 1
pk_guide = 1
pk_tgt_transc = 1

with open(wdir+'/gene_data') as f:
    for linecnt, line in tqdm(enumerate(f)):
        if linecnt % 1000 == 0:
            print("Processing line #%d..."%(linecnt+1))
        if line[0] == '#':
            continue
        geneid, transcid, genename, chrname, desc, cds_start, cds_end = line.strip().split('\t')
        if chrom_table:
            chrname = get_chrom_name(chrname, chrom_table=chrom_table)
        else:
            chrname = get_chrom_name(chrname, chrom_rules=chrom_rules)
        if chrname in err_chroms:
            continue
        cds_start = int(cds_start)-1
        cds_end = int(cds_end)
        try:
            cds_seq = str(chrom_data[chrname][cds_start:cds_end]).upper()
        except KeyError:
            if not chrname in err_chroms:
                print("Couldn't find %s in fasta files! Skipping..."%chrname)
                err_chroms.append(chrname)
            continue

        info_gene = info_genes[geneid]
        if info_gene.pk == -1:
            info_gene.pk = pk_gene
            info_gene.geneid = geneid
            info_gene.name = genename
            info_gene.desc = desc
            pk_gene += 1

        info_transc = info_transcripts[transcid]
        if info_transc.pk == -1:
            info_transc.pk = pk_transc
            info_transc.transcid = transcid
            info_transc.chrname = chrname
            info_transc.gene = info_gene.pk
            info_gene.transcripts.append(transcid) # Assuming no duplication
            pk_transc += 1

        if not (cds_start, cds_end) in info_transc.cds_positions:
            info_transc.cds_positions.append( (cds_start, cds_end) ) # Internally Used

        seqs = []
        for pos in findall_everything(r'[AGTC]{5}[G]{2}', cds_seq):
            i = cds_start + pos - 16
            found_seq = str(chrom_data[chrname][i:i+23]).upper()
            seqs.append( (found_seq, i, '+') )

        for pos in findall_everything(r'[C]{2}[AGTC]{5}', cds_seq):
            i = cds_start + pos
            found_seq = str(chrom_data[chrname][i:i+23].reverse_complement()).upper()
            seqs.append( (found_seq, i, '-') )

        for seq, tgt_pos, direction in seqs:
            if not 'N' in seq:
                pam = seq[-3:]
                seq = seq[:-3]
                targets.append(seq + "NNN 2")

                info_guide = info_guides[seq]
                if info_guide.pk == -1:
                    info_guide.pk = pk_guide
                    info_guide.seq = seq
                    info_guide.gc = GC(seq)
                    pk_guide += 1

                tgt_id = ','.join( [chrname,str(tgt_pos),direction] )
                info_tgt = info_targets[tgt_id]
                if info_tgt.pk == -1:
                    info_tgt.pk = pk_tgt
                    info_tgt.guide = info_guide.pk
                    info_tgt.pam = pam
                    info_tgt.chrname = chrname
                    info_tgt.pos = tgt_pos
                    info_tgt.direction = direction
                    if direction == "+":
                        bp = tgt_pos + len(seq) - 3 # right side of bp
                    else:
                        bp = tgt_pos + 5 # left side of bp

                    if bp-30 < 0 or bp+30 > len(chrom_data[chrname]):
                        info_tgt.oof_score = -1
                    else:
                        info_tgt.oof_score = calc_mich_score(str(chrom_data[chrname][bp-30:bp+30]))[2]
                    pk_tgt += 1

                if not tgt_id in info_transc.targets:
                    info_transc.targets.append(tgt_id)

                if not info_transc.pk in info_tgt.transcripts:
                    info_tgt.transcripts.append(info_transc.pk)


print('Calculating CDS position and coverage of each target...')
for i in tqdm(info_genes):
    for j in info_genes[i].transcripts:
        info_transc = info_transcripts[j]
        tot_cds_length = 0
        for cds_start, cds_end in info_transc.cds_positions:
            tot_cds_length += cds_end - cds_start

        for k in info_transc.targets:
            info_tgt = info_targets[k]
            info_tgt_transc = info_target_transcripts[str(info_tgt.pk)+','+str(info_transc.pk)]
            info_tgt_transc.pk = pk_tgt_transc
            info_tgt_transc.target = info_tgt.pk
            info_tgt_transc.transcript = info_transc.pk
            pk_tgt_transc += 1
            cds_length = 0
            info_transc.cds_positions.sort(key=lambda e: e[0])
            for cds_start, cds_end in info_transc.cds_positions:
                if info_tgt.direction == '+':
                    bp = info_tgt.pos + 20 - 3 # right side of bp
                else:
                    bp = info_tgt.pos + 6 # right side of bp

                if cds_start < bp < cds_end:
                    cds_length += bp - cds_start
                    info_tgt_transc.cds_pos = (cds_length-1) / float(tot_cds_length) * 100.0 # calculated position cannot be the last seq position, thus cds_pos cannot be 100%.
                    info_tgt.coverage += 1
                    break
                else:
                    cds_length += cds_end - cds_start


print('Removing duplicate targets...')
targets = list(set(targets))

print('Writing to file...')
targets_per_file = int(math.ceil(len(targets) / float(TARGET_NUM)))
for cnt, i in enumerate(range(0, len(targets), targets_per_file), 1):
    with open(wdir+'/targets_'+str(cnt)+'.txt', 'w') as fo:
        fo.write(cdir+'\n')
        fo.write('N'*21 + 'RG\n')
        fo.write('\n'.join(targets[i:i+targets_per_file]))

with open(wdir+'/info_genes.txt', 'w') as fo:
    for pk in info_genes:
        gene = info_genes[pk]
        fo.write('\t'.join( [str(gene.pk), gene.geneid, gene.name, gene.desc] ) + "\n" )
with open(wdir+'/info_transcripts.txt', 'w') as fo:
    for pk in info_transcripts:
        transc = info_transcripts[pk]
        fo.write('\t'.join( [str(transc.pk), transc.transcid, transc.chrname, str(transc.gene)] ) + '\t' + '\t'.join( [str(s)+','+str(e) for s,e in transc.cds_positions] ) + "\n" )
with open(wdir+'/info_targets.txt', 'w') as fo:
    for pk in info_targets:
        tgt = info_targets[pk]
        fo.write('\t'.join( [str(tgt.pk), str(tgt.guide), tgt.pam, tgt.chrname, str(tgt.pos), tgt.direction, str(tgt.oof_score), str(tgt.coverage)] ) + '\t' + '\t'.join([str(i) for i in tgt.transcripts]) + '\n')
with open(wdir+'/info_guides.txt', 'w') as fo:
    for pk in info_guides:
        guide = info_guides[pk]
        fo.write('\t'.join( [str(guide.pk), guide.seq, str(guide.gc)] ) + '\n')
with open(wdir+'/info_target_transcripts.txt', 'w') as fo:
    for pk in info_target_transcripts:
        tgt_transc = info_target_transcripts[pk]
        fo.write('\t'.join( [str(tgt_transc.pk), str(tgt_transc.target), str(tgt_transc.transcript), str(tgt_transc.cds_pos)] ) + '\n')

