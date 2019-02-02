#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GO.snakemake.py
@time: 2019/01/31
'''

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]

method=config['diff_para']['methods']
recompare=config['diff_para']['repcompare']
group=config['diff_para']['group']

rule all:
    input:
        expand('{outdir}/GO_Hypergeometric/tmp_file/{recompare}.DEseq2_Method.glist',outdir=outdir,recompare=recompare),
        expand('{outdir}/GO_Hypergeometric/tmp_file/{recompare}.DEseq2_Method.glist2', outdir=outdir, recompare=recompare),
        expand('{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_BP_EnrichResult.sorted.padjust.xls', outdir=outdir,
               recompare=recompare),
        expand('{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_CC_EnrichResult.sorted.padjust.xls', outdir=outdir,
               recompare=recompare),
        expand('{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_MF_EnrichResult.sorted.padjust.xls', outdir=outdir,
               recompare=recompare)

rule goenrichment:
    input:
        '{outdir}/GeneDiffExp_Allin/DEseq2/{recompare}.DEseq2_Method.GeneDiffExpFilter.xls'
    output:
        '{outdir}/GO_Hypergeometric/tmp_file/{recompare}.DEseq2_Method.glist',
        '{outdir}/GO_Hypergeometric/tmp_file/{recompare}.DEseq2_Method.glist2',
        '{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_BP_EnrichResult.sorted.padjust.xls',
        '{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_CC_EnrichResult.sorted.padjust.xls',
        '{outdir}/GO_Hypergeometric/GO/{recompare}.DEseq2_Method_MF_EnrichResult.sorted.padjust.xls'
    params:
        goclass="/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/go.class",
        annot="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/GO/database/hg19.annot",
        outdir="{outdir}/GO_Hypergeometric/GO",
        tmpdir="{outdir}/GO_Hypergeometric/tmp_file",
        species="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/GO/GOEnrichment/species",
        obo="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/GO/GOEnrichment/gene_ontology.1_2.obo",
        sb="/ldfssz1/ST_BIGDATA/USER/yueyao/bin/GO/GOEnrichment",
        py="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/envs/snakemake/bin/python"
    shell:
        "awk '{{print $1\"\\t\"$7}}' {input} >{params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist;" \
        "awk '{{print $1}}' {input} >{params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist2;" \
        "sed -i '1d' {params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist2;" \
        "sed -i '1d' {params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist;" \
        "{params.py} {params.sb}/drawGO.py --goclass {params.goclass} --annot {params.annot} --gene {params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist2 --outprefix {params.outdir}/{wildcards.recompare}.DEseq2_Method;" \
        "{params.py} {params.sb}/GOTermEnrichment.py --oboFile {params.obo} --selectionSetFiles {params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist2 --gene2GoFile {params.annot} --outdir {params.outdir}; " \
        "{params.py} {params.sb}/topGO.py --godir {params.outdir} --glist {params.tmpdir}/{wildcards.recompare}.DEseq2_Method.glist --outdir {params.outdir} --prefix {params.species}"