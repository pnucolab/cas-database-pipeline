#!/bin/bash
nohup pypy 2.generate_targets.py athaliana_eg_gene arabidopsis_TAIR10 > stdout_gentarget_at.txt &

nohup pypy 2.generate_targets.py celegans_gene_ensembl C.elegans_ce10 > stdout_gentarget_ce.txt &

nohup pypy 2.generate_targets.py dmelanogaster_gene_ensembl Drosophila_dm3 > stdout_gentarget_dm.txt &

nohup pypy 2.generate_targets.py drerio_gene_ensembl zebrafish_danRer7 > stdout_gentarget_dr.txt &

nohup pypy 2.generate_targets.py hsapiens_gene_ensembl human_hg38 > stdout_gentarget_hs.txt &

nohup pypy 2.generate_targets.py mmusculus_gene_ensembl mouse_mm10 > stdout_gentarget_mm.txt &

nohup pypy 2.generate_targets.py rnorvegicus_gene_ensembl Rat_rn5 > stdout_gentarget_rn.txt &

nohup pypy 2.generate_targets.py slycopersicum_eg_gene tomato_SL2.4 > stdout_gentarget_sl.txt &

nohup pypy 2.generate_targets.py xtropicalis_gene_ensembl xenopus42 > stdout_gentarget_xt.txt &
