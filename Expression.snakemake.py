#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Expression.snakemake.py
@time: 2019/01/09
'''

include: 'Prediction.snakemake.py'

SAMPLES = config["samples"]
soft=config["software"]
sb=config["software"]["scriptbin"]
outdir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3"

rule all:
    input:
        expand('{outdir}/Expression/ciri_multisamples_readscount.txt',outdir=outdir),
        expand('{outdir}/Expression/circRNA_expression.xls',outdir=outdir),
        expand('{outdir}/Expression/multisamples_strand_count.txt',outdir=outdir),
        expand('{outdir}/Expression/Boxplot.png',outdir=outdir)

rule expression:
    input:
        ciri=expand('{outdir}/Prediction/{sample}/{sample}_ciri_output.xls',outdir=outdir,sample=SAMPLES),
        findciri=expand('{outdir}/Prediction/{sample}/{sample}_findcirc_output.xls',outdir=outdir,sample=SAMPLES),
        bowtie2log=expand('{outdir}/Prediction/{sample}/{sample}_bowtie2.log',outdir=outdir,sample=SAMPLES)
    output:
        '{outdir}/Expression/ciri_multisamples_readscount.txt',
        '{outdir}/Expression/findcirc_multisamples_readscount.txt',
        '{outdir}/Expression/circRNA_expression.xls',
        '{outdir}/Expression/Boxplot.png',
        '{outdir}/Expression/circID_linktable.txt',
        '{outdir}/Expression/multisamples_strand_count.txt'
    params:
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Expression/"
    run:
        cirilist=",".join(input.ciri)
        findcirilist=",".join(input.findciri)
        bowtie2loglist=",".join(input.bowtie2log)
        od=os.path.dirname(output[0])
        shell("perl {params.sb}/multisamples_readscount.pl -ciri {cirilist} -findcirc {findcirilist} -outdir {od}; \
            perl {params.sb}/circID_linking.pl -ciri {od}/ciri_multisamples_readscount.txt -findcirc {od}/findcirc_multisamples_readscount.txt -outdir {od}; \
            perl {params.sb}/calculate_expression.pl -ciri {od}/ciri_multisamples_readscount.txt -findcirc {od}/findcirc_multisamples_readscount.txt -link {output[4]} -bowtielog {bowtie2loglist} -outdir {od}; \
            perl {params.sb}/drawbox_density-plot.pl -i  {od} -o {od}")

