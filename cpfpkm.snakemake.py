#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: cpfpkm.snakemake.py
@time: 2018/12/28
'''

sample=['HBRR1','HBRR2','HBRR3','UHRR1','UHRR2','UHRR3']
outdir='/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake_md/test/snakemake_md-example/result/BGI_result/Quantify/GeneExpression/GeneExpression'
indir='/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake_md/test/snakemake_md-example/result'

rule all:
    input:
        expand('{outdir}/{sample}.gene.fpkm.xls',outdir=outdir,sample=sample)

# rule cp:
#     input:
#         expand('{indir}/GeneExp_RSEM/{sample}/{sample}.gene.fpkm.xls',indir=indir,sample=sample)
#     output:
#         '{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls'
#     run:
#         id=indir
#         od=outdir
#         shell('cp {id}/GeneExp_RSEM/{wildcards.sample}/{wildcards.sample}.gene.fpkm.xls {od}/GeneExp/{wildcards.sample}/{wildcards.sample}.gene.fpkm.xls')
#
# rule cp2:
#     input:
#         expand('{indir}/{sample}.gene.fpkm.xls',indir=indir,sample=sample)
#     output:
#         '{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls'
#     run:
#         id=indir
#         od=outdir
#         shell('cp {id}/{wildcards.sample}.gene.fpkm.xls {od}/GeneExp/{wildcards.sample}/{wildcards.sample}.gene.fpkm.xls')

rule cp3:
    input:
        expand('{indir}/GeneExp/{sample}/{sample}.gene.fpkm.xls',indir=indir,sample=sample)
    output:
        '{outdir}/{sample}.gene.fpkm.xls'
    run:
        id=indir
        od=outdir
        shell('cp {id}/GeneExp/{wildcards.sample}/{wildcards.sample}.gene.fpkm.xls {od}/{wildcards.sample}.gene.fpkm.xls')
