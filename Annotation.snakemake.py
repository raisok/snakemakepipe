#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Annotation.snakemake.py
@time: 2019/01/09
'''

SAMPLES = config["samples"]
soft=config["software"]
sb=config["software"]["scriptbin"]
outdir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3"

rule all:
    input:
        expand('{outdir}/Annotation/{sample}_integrated_prediction.xls',outdir=outdir,sample=SAMPLES),
        expand('{outdir}/Annotation/barplot/circBase_num.xls',outdir=outdir),
        expand('{outdir}/Annotation/barplot/circRNA_typenum.xls',outdir=outdir),
        expand('{outdir}/Annotation/Venn/{sample}_integrated_prediction.venn.png',outdir=outdir,sample=SAMPLES),
        expand("{outdir}/Annotation/{sample}_integrated_prediction_desc.xls", outdir=outdir, sample=SAMPLES)

# circRNA circBase_ID, Strand, Type
rule annotation_step01:
    input:
        expand("{outdir}/Expression/multisamples_strand_count.txt",outdir=outdir)
    output:
        expand("{outdir}/Annotation/{sample}_integrated_prediction.xls",outdir=outdir,sample=SAMPLES),
        expand('{outdir}/Annotation/barplot/circBase_num.xls',outdir=outdir),
        expand('{outdir}/Annotation/barplot/circRNA_typenum.xls',outdir=outdir)
    params:
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Annotation/",
        base="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/database/circBase/hsa_hg19_circRNA.txt",
        gtf="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation/hg19_filter.psl.gtf",
        od=expand("{outdir}/Annotation/",outdir=outdir)
    shell:
        "perl {params.sb}/circRNA_type.pl {params.base} {input} {params.gtf} {params.od}; \
         perl {params.sb}/barplot.circBaseNumber.pl -indir {params.od} -outdir {params.od}/barplot; \
         perl {params.sb}/barplot.TypeNumber.pl -indir {params.od} -outdir {params.od}/barplot"

# circRNA prediction Venn Diagram
rule annotation_step02:
    input:
        '{outdir}/Annotation/{sample}_integrated_prediction.xls'
    output:
        '{outdir}/Annotation/Venn/{sample}_integrated_prediction.venn.png'
    params:
        od='{outdir}/Annotation/Venn/',
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Annotation/"
    shell:
        "perl {params.sb}/venny.pl -infile {input} -name CIRI,find_circ -header -outdir {params.od} -imgname {wildcards.sample}_integrated_prediction.venn"


# circRNA Annotation
rule annotation_step03:
    input:
        expand("{outdir}/Annotation/{sample}_integrated_prediction.xls", outdir=outdir, sample=SAMPLES)
    output:
        expand("{outdir}/Annotation/{sample}_integrated_prediction_desc.xls",outdir=outdir,sample=SAMPLES)
    params:
        gene2mark='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg76/refMrna.fa.gene2mark',
        gene2symbol='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/human.gene2symbol.txt',
        go='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19',
        kegg='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19.ko',
        nr='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19.nr.desc',
        db='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_2016a/prepare/NovelTr_db_version.txt',
        spe='an',
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Annotation",
        od=expand('{outdir}/Annotation',outdir=outdir)
    shell:
        "perl {params.sb}/genDesc.pl {params.gene2mark} {params.gene2symbol} {params.go} {params.kegg} {params.nr} {params.db} {params.spe} {params.od}/Description.xls; \
         perl {params.sb}/circRNA_addDesc.pl {params.od}/Description.xls {params.od} {params.od}"