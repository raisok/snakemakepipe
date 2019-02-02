#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: CircosFig.snakemake.py
@time: 2019/01/09
'''

outdir='/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3'
SAMPLES = config["samples"]

rule all:
    input:
        expand('{outdir}/CircosFig/{sample}/chr_basenum.txt',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/CircosFig/{sample}/karyotype.txt',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/CircosFig/{sample}/{sample}.circRNAnum.txt',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/CircosFig/{sample}/{sample}.circRNAnum.circos.txt',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/CircosFig/{sample}/{sample}.circos.svg',outdir=outdir,sample=SAMPLES),

rule circos_step01:
    input:
        '/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/prepare/GenomeBowtie2Index/chrALL.fa'
    output:
        '{outdir}/CircosFig/{sample}/chr_basenum.txt'
    params:
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Circos"
    shell:
        "perl {params.sb}/seq_basenum.pl {input} {output}"

rule circos_step02:
    input:
        rules.circos_step01.output

    output:
        '{outdir}/CircosFig/{sample}/karyotype.txt'
    params:
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Circos/"
    shell:
        "perl {params.sb}/circos_karyotype.pl {input} {output}"

rule circos_step03:
    input:
        rules.circos_step01.output,
        '{outdir}/Annotation/{sample}_integrated_prediction.xls'
    output:
        '{outdir}/CircosFig/{sample}/{sample}.circRNAnum.txt'
    shell:
        "grep -v '^Final_ID' {input[1]} | awk '{{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4}}' >{output}"

rule circos_step04:
    input:
        rules.circos_step03.output,
        rules.circos_step01.output
    output:
        '{outdir}/CircosFig/{sample}/{sample}.circRNAnum.circos.txt'
    params:
        sb = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Circos/",
        od="{outdir}/CircosFig/{sample}"
    shell:
        "perl {params.sb}/circRNA_num.pl {input[0]} {input[1]} {params.od} 1000000; \
        cat {params.od}/{wildcards.sample}.circRNAnum.*.xls >{output}"

rule circos_step05:
    input:
        rules.circos_step02.output,
        rules.circos_step04.output
    output:
        '{outdir}/CircosFig/{sample}/{sample}.circos.svg'
    params:
        sb = "/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Circos/",
        od="{outdir}/CircosFig/{sample}"
    shell:
        "perl {params.sb}/draw_circRNA_circos.pl --kary {input[0]} --unit 1000000 --circ {input[1]}  --label 10p --tick 5p --output {params.od}"