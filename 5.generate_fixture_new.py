import os, sys, django
sys.path.append('/data/django/rgenome_net')
os.environ['DJANGO_SETTINGS_MODULE'] = 'rgenome_net.settings'
from django.conf import settings
django.setup()

from cas_database.models import *

from sys import argv
import json
import codecs

p = argv[1].rstrip('/')

org_id = Organism.objects.all().count() + 1
org_str = ""
while org_str == "":
    org_str = raw_input("Enter Organism Name (e.g. Homo sapiens (GRCh38/hg38) - Human):\n")

print("")
for ot in OrganismType.objects.all():
    print("{0}. {1}".format(ot.pk, ot.name))

org_type = 0
while org_type == 0:
    try:
        org_type = int(raw_input("Enter Organism Type:\n"))
    except:
        org_type = 0

def load_file(filename, info_dic):
    with open(p+'/'+filename) as f:
        for cnt, line in enumerate(f, start=1):
            entries = line.strip().split('\t')
            info_dic[int(entries[0])] = entries[1:]

class fixture:
    def __init__(self, t, p):
        if not os.path.isdir(p+"/fixtures"):
            os.mkdir(p+"/fixtures")
        self.fixture_item_count = 0
        self.fixture_count = 1
        self.fixture_file = codecs.open(p+"/fixtures/fixture_%s_00001.json"%t, 'w', "utf-8")
        self.fixture_file.write('[')
        self.t = t

    def __enter__(self):
        return self

    def append(self, s, noend=False):
        if not noend and self.fixture_item_count >= 1000000:
            self.fixture_file.write(']')
            self.fixture_file.close()
            self.fixture_count += 1
            self.fixture_file = open(p+"/fixtures/fixture_%s_%05d.json"%(self.t, self.fixture_count), 'w')
            self.fixture_file.write('[')
            self.fixture_item_count = 0
        else:
            if self.fixture_item_count != 0:
                self.fixture_file.write(',\n')
        self.fixture_file.write(s)
        self.fixture_item_count += 1

    def close(self):
        self.fixture_file.write(']')
        self.fixture_file.close()

    def __exit__(self, t, v, tr):
        try:
            self.close()
        except:
            pass

print("Loading off-target information into memory...")
info_offtargets = {}
load_file('info_cpp_offtargets.txt', info_offtargets)

with fixture("org", p) as f:
    print("Generating fixtures (Organism)...")
    f.append( json.dumps( { "model": "cas_database.organism",
                            "pk": org_id,
                            "fields": { "name": org_str,
                                        "assemblyid": "",
                                        "org_type": org_type } } ) )

cnt_gene = Gene.objects.all().count()
with fixture("gene", p) as f:
    print("Generating fixtures (Genes)...")
    with open(p+"/info_genes.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, gene = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.gene",
                                    "pk": cnt_gene + pk,
                                    "fields": { "organism": org_id,
                                                "symbol": gene[1],
                                                "geneid": gene[0],
                                                "name": gene[2] } } ) )

info_chroms = {}
pk_chrom = Chromosome.objects.all().count() + 1
cnt_transcript = Transcript.objects.all().count()
cnt_cdsrange = CDSrange.objects.all().count()
with fixture("ctc", p) as f:
    print("Generating fixtures (Chromosome, Transcripts, and CDSrange)...")
    pk_range = CDSrange.objects.all().count() + 1
    with open(p+"/info_transcripts.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, transc = int(entries[0]), entries[1:]
            if not transc[1] in info_chroms:
                f.append( json.dumps( { "model": "cas_database.chromosome",
                                        "pk": pk_chrom,
                                        "fields": { "organism": org_id,
                                                    "name": transc[1],
                                                    "ncbiname": "" } } ) )
                info_chroms[transc[1]] = pk_chrom
                pk_chrom += 1

            f.append( json.dumps( { "model": "cas_database.transcript",
                                    "pk": cnt_transcript + pk,
                                    "fields": { "organism": org_id,
                                                "transcriptid": transc[0],
                                                "chrom": info_chroms[transc[1]],
                                                "gene": cnt_gene + int(transc[2]),
                                                "max_coverage": 0 } } ), True )
            for r in transc[3:]:
                a, b = r.split(",")[:2]
                f.append( json.dumps( { "model": "cas_database.cdsrange",
                                        "pk": cnt_cdsrange + pk_range,
                                        "fields": { "transcript": cnt_transcript + pk,
                                                    "start": int(a),
                                                    "end": int(b) } } ), True )
                pk_range += 1

print("Loading all sequences in database...")
cnt_sequence = Sequence.objects.all().count()
seqs = {i.seq: i.id for i in Sequence.objects.all()}

cnt_guide = Guide.objects.all().count()
with fixture("gs", p) as f:
    print("Generating fixtures (Guides and Sequences)...")
    with open(p+"/info_guides.txt") as fi, open(p+"/targets_notfound.txt", "w") as fo:
        fo.write('human_hg38\n'+'N'*21+'RG\n')
        for line in fi:
            entries = line.strip().split('\t')
            pk, seq = int(entries[0]), entries[1:]
            try:
                off = info_offtargets[pk]
            except:
                fo.write(seq[0] + "NNN 2\n")
            try:
                seq_key = seqs[seq[0]]
            except:
	        cnt_sequence += 1
                seq_key = cnt_sequence
                f.append( json.dumps( { "model": "cas_database.sequence",
                                        "pk": seq_key,
                                        "fields": { "seq": seq[0],
                                                    "gc": float(seq[1]) } } ) )

            f.append( json.dumps( { "model": "cas_database.guide",
                                    "pk": cnt_guide + pk,
                                    "fields": { "organism": org_id,
                                                "seq": seq_key,
                                                "off_0": int(off[0]),
                                                "off_1": int(off[1]),
                                                "off_2": int(off[2]) } } ), True )

cnt_target = Target.objects.all().count()
with fixture("tgt", p) as f:
    print("Generating fixtures (Targets)...")
    with open(p+"/info_targets.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, tgt = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.target",
                                    "pk": cnt_target + pk,
                                    "fields": { "organism": org_id,
                                                "guide": cnt_guide + int(tgt[0]),
                                                "pam": tgt[1],
                                                "chrom": info_chroms[tgt[2]],
                                                "loci": int(tgt[3]),
                                                "rev": True if tgt[4] == "-" else False,
                                                "oof_score": float(tgt[5]),
                                                "cov": int(tgt[6]) } } ) )

cnt_target_transcript = Target_Transcript.objects.all().count()
with fixture("tt", p) as f:
    print("Generating fixtures (Target-Transcripts)...")
    with open(p+"/info_target_transcripts.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, tgt_transc = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.target_transcript",
                                    "pk": cnt_target_transcript + pk,
                                    "fields": { "target": cnt_target + int(tgt_transc[0]),
                                                "transcript": cnt_transcript + int(tgt_transc[1]),
                                                "cds_percentage": float(tgt_transc[2])} } ) )

cnt_offtarget = Offtarget.objects.all().count()
with fixture("ot", p) as f:
    print("Generating fixtures (Off-Targets)...")
    with open(p+"/info_cpp_offtargets_ref.txt") as fi:
        for line in fi:
            entries = line.strip().split('\t')
            pk, off = int(entries[0]), entries[1:]
            f.append( json.dumps( { "model": "cas_database.offtarget",
                                    "pk": cnt_offtarget + pk,
                                    "fields": { "guide": cnt_guide + int(off[0]),
                                                "off_count": int(off[1]),
                                                "chrom": info_chroms[off[2]],
                                                "loci": int(off[3]),
                                                "seq": off[4],
                                                "rev": True if tgt[4] == "-" else False,
                                                "transcript": None if off[6] == "-1" else cnt_transcript + int(off[6]) } } ) )
