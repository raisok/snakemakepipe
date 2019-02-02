#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: DEG.circosRNA.snakemake.py
@time: 2019/01/10
'''

SAMPLES = config["samples"]
soft=config["software"]
outdir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3"
sb=config["software"]["scriptbin"]

method=config['diff_para']['methods']
repcompare=config['diff_para']['repcompare']
norepcompare=config['diff_para']['norepcompare']
group=config['diff_para']['group']

rule all:
    input:
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.DiffExp.xls',outdir=outdir,norepcompare=norepcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.DiffExpFilter.xls',outdir=outdir,norepcompare=norepcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.Scatter-plot.png',outdir=outdir,norepcompare=norepcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.Volcano-plot.png',outdir=outdir,norepcompare=norepcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.pheatmap-plot.png',outdir=outdir,norepcompare=norepcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExp.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExpFilter.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.Scatter-plot.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.Volcano-plot.png',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.pheatmap-plot.png',outdir=outdir,repcompare=repcompare)

rule pre:
    input:
        expand('{outdir}/Expression/circRNA_expression.xls',outdir=outdir)
    output:
        difflist='{outdir}/GeneDiffExp_Allin/temp/diffCompare',
        grouplist='{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist='{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf='{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    params:
        sample=SAMPLES
    run:
        redir = os.path.dirname(output.difflist)
        os.makedirs(redir,exist_ok=True)
        group_str=""
        with open(redir+'/GroupList.txt','w') as gl:
            for j in group.keys():
                gl.write(j+"\t"+group[j]+"\n")
                group_str+=j+":"+group[j]+" "
        group_str=group_str.rstrip()
        repcompare_str=','.join(repcompare).replace('-VS-','&')
        norepcompare_str=','.join(norepcompare).replace('-VS-','&')
        with open(redir+'/conf.lst','w') as conf:
            for i in method.keys():
                if i == 'DEGseq':
                    conf.write('GeneDiff_Method = ' + i + '\n')
                    conf.write(i + '_VS = ' + norepcompare_str + '\n')
                    conf.write(i + '_Filter = ' + method[i] + '\n\n')
                elif i == 'DEseq2':
                    conf.write('GeneDiff_Method = '+i+'\n')
                    conf.write(i+'_Group = '+group_str+'\n')
                    conf.write(i+ '_VS = '+repcompare_str+'\n')
                    conf.write(i + '_Filter = ' + method[i] + '\n\n')
        with open (redir+'/CompareList.txt','w') as cl:
            for clt in repcompare:
                arr = clt.split('-VS-')
                cl.write(group[arr[0]]+"\t"+group[arr[1]]+"\n")
        with open(redir+'/diffCompare','w') as dc:
            for mm in norepcompare:
                mm=mm.replace("-VS-",'\t')
                dc.write(mm+'\n')


rule degseq:
    input:
        exp='{outdir}/Expression/circRNA_expression.xls',
        difflist='{outdir}/GeneDiffExp_Allin/temp/diffCompare',
        conf='{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.DiffExp.xls',
        '{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.DiffExpFilter.xls',
        '{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.Scatter-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.Volcano-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.pheatmap-plot.png'
    params:
        script = '/ldfssz1/ST_BIGDATA/USER/yueyao/bin/RNA_Script/CirRNA/DEG',
        para = method['DEGseq']
    run:
        outd = os.path.dirname(output[0])
        shell("perl {params.script}/DEGseq.pl -count {input.exp} -diff {input.difflist} -method MARS -threshold 5 -pValue 1e-3 -zScore 4 {params.para} -outdir {outd}; \
              perl {params.script}/drawScatter-plot.pl -indir {outd} -exp1col 2 -exp2col 3 -udcol 6 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/drawVolcano-plot.pl -Indir {outd} -log2col 4 -signcol 5 -udcol 6 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/draw_pheatmap.pl -diffdir {outd} -explist {input.exp} -config {input.conf} -outdir {outd} -udcol 6")



rule deseq2:
    input:
        exp = '{outdir}/Expression/circRNA_expression.xls4',
        grouplist = '{outdir}/GeneDiffExp_Allin/temp/GroupList.txt',
        comparelist = '{outdir}/GeneDiffExp_Allin/temp/CompareList.txt',
        conf = '{outdir}/GeneDiffExp_Allin/temp/conf.lst'
    output:
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExp.xls',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExpFilter.xls',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.Scatter-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.Volcano-plot.png',
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.pheatmap-plot.png'
    params:
        script='/ldfssz1/ST_BIGDATA/USER/yueyao/bin/RNA_Script/CirRNA/DEG',
        para=method['DEseq2']
    run:
        outd=os.path.dirname(output[0])
        shell("perl {params.script}/DEseq2.circRNA.pl -count {input.exp} -diff {input.comparelist} -group {input.grouplist} {params.para} -outdir {outd}; \
              perl {params.script}/drawScatter-plot.pl -indir {outd} -exp1col 2 -exp2col 3 -udcol 6 -conf {input.conf} -outdir {outd}; \
              perl {params.script}/drawVolcano-plot.pl -Indir {outd} -log2col 4 -signcol 5 -udcol 6 -xlab \"log2(fold change)\" -ylab=\"-log10(Padj)\" -conf {input.conf} -outdir {outd}; \
              perl {params.script}/draw_pheatmap.pl -diffdir {outd} -explist {input.exp} -config {input.conf} -outdir {outd} -udcol 6")
