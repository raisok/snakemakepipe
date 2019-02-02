#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os


'''
@author: yueyao
@file: GeneDiff.snakemake_md.py
@time: 2018/12/27
'''

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]

method=config['diff_para']['methods']
compare=config['diff_para']['repcompare']
group=config['diff_para']['group']

# rule all:
#     input:
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExp.xls',outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.GeneDiffExp.xls',outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.GeneDiffExp.xls',outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExpFilter.xls', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.GeneDiffExpFilter.xls', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.GeneDiffExpFilter.xls', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.MA-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.MA-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.MA-plot.png', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.Scatter-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.Scatter-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.Scatter-plot.png', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.Volcano-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.Volcano-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.Volcano-plot.png', outdir=outdir, compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.pheatmap-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.pheatmap-plot.png', outdir=outdir,compare=compare),
#         expand('{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.pheatmap-plot.png', outdir=outdir, compare=compare)

rule pre:
    input:
        expand('{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls',outdir=outdir,sample=SAMPLES)
    output:
        samplelist='{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist='{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist='{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf='{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    params:
        sample=SAMPLES
    run:
        redir = os.path.dirname(output.samplelist)
        os.makedirs(redir,exist_ok=True)
        group_str=""
        with open(redir+'/GroupList.txt','w') as gl:
            for j in group.keys():
                gl.write(j+"\t"+group[j]+"\n")
                group_str+=j+":"+group[j]+" "
        group_str=group_str.rstrip()
        compare_str=','.join(compare).replace('-VS-','&')
        with open(redir+'/conf.lst','w') as conf:
            for i in method.keys():
                conf.write('GeneDiff_Method = '+i+'\n')
                conf.write(i+'_Group = '+group_str+'\n')
                conf.write(i+ '_VS = '+compare_str+'\n')
                conf.write(i + '_Filter = ' + method[i] + '\n\n')
        with open(redir+'/SampleList.txt','w') as sl:
            for fpkm in input:
                name=os.path.basename(fpkm).split('.')[0]
                sl.write(name+'\t'+fpkm+"\n")
        with open (redir+'/CompareList.txt','w') as cl:
            for clt in compare:
                arr = clt.split('-VS-')
                cl.write(group[arr[0]]+"\t"+group[arr[1]]+"\n")

rule deseq2:
    input:
        samplelist = '{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist = '{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExp.xls',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExpFilter.xls',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.MA-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.Scatter-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.Volcano-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.pheatmap-plot.png'

    params:
        script=sb[0],
        para=method['DEseq2']
    run:
        outd=os.path.dirname(output[0])
        shell("perl {params.script}/GeneDiffExp/DEseq2.pl -list {input.samplelist} -diff {input.comparelist} -group {input.grouplist} {params.para} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawMA-plot.pl -indir {outd}  -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawScatter-plot.pl -indir {outd} -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawVolcano-plot.pl -Indir {outd} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/draw_pheatmap.pl -diffdir {outd} -explist {input.samplelist} -config {input.conf} -colgeneid 1 -colexp 5 -outdir {outd} -udcol 7")


rule noiseq:
    input:
        samplelist = '{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist = '{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.GeneDiffExp.xls',
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.GeneDiffExpFilter.xls',
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.MA-plot.png',
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.Scatter-plot.png',
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.Volcano-plot.png',
        '{outdir}/GeneDiffExp_Allin/NOIseq/{compare}.NOIseq_Method.pheatmap-plot.png'
    params:
        script = sb[0],
        para = method['NOIseq']
    run:
        outd = os.path.dirname(output[0])
        shell("perl {params.script}/GeneDiffExp/NOIseq.pl -list {input.samplelist} -diff {input.comparelist} -group {input.grouplist} {params.para} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawMA-plot.pl -indir {outd}  -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawScatter-plot.pl -indir {outd} -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawVolcano-plot.pl -Indir {outd} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/draw_pheatmap.pl -diffdir {outd} -explist {input.samplelist} -config {input.conf} -colgeneid 1 -colexp 5 -outdir {outd} -udcol 7")

rule ebseq:
    input:
        samplelist = '{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist = '{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.GeneDiffExp.xls',
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.GeneDiffExpFilter.xls',
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.MA-plot.png',
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.Scatter-plot.png',
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.Volcano-plot.png',
        '{outdir}/GeneDiffExp_Allin/EBseq/{compare}.EBseq_Method.pheatmap-plot.png'
    params:
        script = sb[0],
        para = method['EBseq']
    run:
        outd = os.path.dirname(output[0])
        shell("perl {params.script}/GeneDiffExp/EBseq.pl -list {input.samplelist} -diff {input.comparelist} -group {input.grouplist} {params.para} -outdir {outd}; \
                  perl {params.script}/GeneDiffExp/drawMA-plot.pl -indir {outd}  -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
                  perl {params.script}/GeneDiffExp/drawScatter-plot.pl -indir {outd} -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
                  perl {params.script}/GeneDiffExp/drawVolcano-plot.pl -Indir {outd} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
                  perl {params.script}/GeneDiffExp/draw_pheatmap.pl -diffdir {outd} -explist {input.samplelist} -config {input.conf} -colgeneid 1 -colexp 5 -outdir {outd} -udcol 7")

