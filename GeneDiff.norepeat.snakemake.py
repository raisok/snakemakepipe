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
compare=config['diff_para']['compare']
group=config['diff_para']['group']

rule all:
    input:
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExp.xls',outdir=outdir,compare=compare)

rule pre:
    input:
        expand('{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls',outdir=outdir,sample=SAMPLES)
    output:
        samplelist='{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist='{outdir}/GeneDiffExp_Allin/temp/diffCompare',
        comparelist='{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf='{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    params:
        sample=SAMPLES
    run:
        redir = os.path.dirname(output.samplelist)
        os.makedirs(redir,exist_ok=True)
        compare_str=','.join(compare).replace('-VS-','&')
        with open(redir+'/conf.lst','w') as conf:
            for i in method.keys():
                conf.write('GeneDiff_Method = '+i+'\n')
                conf.write(i+ '_VS = '+compare_str+'\n')
                conf.write(i + '_Filter = ' + method[i] + '\n\n')
        with open(redir+'/SampleList.txt','w') as sl:
            for fpkm in input:
                dic_c[name] = fpkm
                name=os.path.basename(fpkm).split('.')[0]
                sl.write(name+'\t'+fpkm+"\n")
        with open (redir+'/CompareList.txt','w') as cl:
            for clt in compare:
                arr = clt.split('-VS-')
                cl.write(arr[0]+'\t'+dic_c[arr[0]]+'\n')
                cl.write(arr[1] + '\t' + dic_c[arr[1]] + '\n')
        with open(redir+'/diffCompare','w') as dc:
            for clt in compare:
                arr = clt.split('-VS-')
                dc.write(arr[0]+'\t'+arr[1]+'\n')


rule possion:
    input:
        samplelist='{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        grouplist = '{outdir}/GeneDiffExp_Allin/temp/diffCompare',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/Possion/{compare}.Possion_Method.GeneDiffExp.xls'
    params:
        script = sb[0],
        para = method['Possion']
    run:
        outd = os.path.dirname(output[0])
        shell("perl {params.script}/GeneDiffExp/Possion.pl -list {input.comparelist} {params.para} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawMA-plot.pl -indir {outd}  -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawScatter-plot.pl -indir {outd} -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawVolcano-plot.pl -Indir {outd} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/draw_pheatmap.pl -diffdir {outd} -explist {input.samplelist} -config {input.conf} -colgeneid 1 -colexp 5 -outdir {outd} -udcol 7")

rule degseq:
    input:
        samplelist='{outdir}/GeneDiffExp_Allin/temp/SampleList.txt',
        difflist = '{outdir}/GeneDiffExp_Allin/temp/diffCompare',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/DEseq2/{compare}.DEseq2_Method.GeneDiffExp.xls'
    params:
        script = sb[0],
        para = method['DEGseq']
    run:
        outd = os.path.dirname(output[0])
        shell("perl {params.script}/GeneDiffExp/DEGseq.pl -list {input.comparelist} -diff {input.difflist} -GeneIDColumn 1 -GeneExpColumn 5 -GeneExpreadscountColumn 4 -GeneLenColumn 3 -method MARS -threshold 5 -pValue 1e-3 -zScore 4 {params.para} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawMA-plot.pl -indir {outd}  -log2col 5 -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawScatter-plot.pl -indir {outd} -exp1col 3 -exp2col 4 -udcol 7 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/drawVolcano-plot.pl -Indir {outd} -log2col 5 -signcol 6 -udcol 7 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/GeneDiffExp/draw_pheatmap.pl -diffdir {outd} -explist {input.samplelist} -config {input.conf} -colgeneid 1 -colexp 5 -outdir {outd} -udcol 7")

