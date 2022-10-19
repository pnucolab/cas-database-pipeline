with open('SLmic1.0_gene_models.gff') as fp:
    data = fp.readlines()
    data = [item.replace('\n','') for item in data]
    for i in range(len(data)):
        data[i] = data[i].split()
cds_dict = {}
for i in range(len(data)):
    #print(data[i][0])
    cds_dict[f'{data[i][0]}'] = []
# print(cds_dict)
for i in range(len(data)):
    if data[i][2] == 'CDS':
        cds_dict[f'{data[i][0]}'].append([data[i][-1][3:13],int(data[i][3])-1,int(data[i][4])])# 파일의 1-based에서 0-based로 바뀐 것을 고려(단, 슬라이싱을 고려하여 end좌표는 1을 빼지 않음)
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

from Bio import SeqIO
# 파일 형식: (sequence), (chr), (location), (+/-), (gene)
cleavage_dict = {}
chromo_list = ['SLmic1.0_chr1','SLmic1.0_chr2','SLmic1.0_chr3','SLmic1.0_chr4','SLmic1.0_chr5','SLmic1.0_chr6','SLmic1.0_chr7','SLmic1.0_chr8','SLmic1.0_chr9','SLmic1.0_chr10','SLmic1.0_chr11','SLmic1.0_chr12']

for i in chromo_list:
    cleavage_dict[i] = []


for seq_record in SeqIO.parse("GCA_012431665.1_SLYMIC_genomic.fa", "fasta"):
    seq = seq_record.seq

    seq=str(seq.upper()) 

    print(seq_record.id)

    # cleavage, start_point, sequence, ngg/ccn list
    if seq_record.id in chromo_list:
        for i in range(21,len(seq)-1):
            if seq[i] == 'G' and seq[i+1]=='G':# ngg-3, ngg-20, sequence
                cleavage_dict[f'{seq_record.id}'].append([i-4,i-21,seq[i-21:i+2],'+'])
        for i in range(len(seq)-24):
            if seq[i] == 'C' and seq[i+1]=='C' and (seq_record.id in chromo_list):# ccn+6, ccn, sequence
                cleavage_dict[f'{seq_record.id}'].append([i+5,i,seq[i:i+23][::1],'-'])


# 파일 형식: (sequence), (chr), (location), (+/-), (gene)
f = open("first.txt", 'w')
f.write('# (sequence), (chr), (location), (+/-), (gene)\n')
for i in cleavage_dict:
    print(i)
    for j in cds_dict[i]:
        for k in cleavage_dict[i]:
            if (j[1] <= k[0]) and (k[0] < j[2]) :
                f.write(f'({k[2]}), ({i}), ({k[1]}), ({k[-1]}), ({j[0]})\n')
f.close()

