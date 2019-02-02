#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: RNAseq.snakemake.py
@time: 2018/12/19
'''

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]
ruledir=config["ruledir"]

filter_step=config["step"]["filter"]
alignment_step=config["step"]["alignment"]
geneexp_step=config["step"]["geneexp"]
genediff_step=config['step']['genediff']
go_step=config['step']['goenrichment']

ifrepeat=config['diff_para']['group']

for i_v,i_k in ifrepeat.items():
    arr=i_k.split(",")
    if len(arr) >1:
        repcompare = config['diff_para']['repcompare']
        rep=1
    if len(arr) ==1:
        repcompare = config['diff_para']['norepcompare']
        rep=0

rule all:
    input:
        expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz', outdir=outdir,sample=SAMPLES),
        expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz', outdir=outdir,sample=SAMPLES),
        # expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.base.png', outdir=outdir,sample=SAMPLES),
        # expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.qual.png', outdir=outdir,sample=SAMPLES),
        # expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.RawReadsClass.png', outdir=outdir,sample=SAMPLES),
        # expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.filter.stat.xls', outdir=outdir,sample=SAMPLES),
        # expand('{outdir}/Filter_SOAPnuke/FilterSummary.xls', outdir=outdir),
        expand('{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam.bai', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GenomeMapping_HISAT/GenomeMappingSummary.xls', outdir=outdir),
        expand('{outdir}/GeneExp/{sample}/{sample}.Map2GeneStat.xls', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GeneExp/{sample}/{sample}.isoforms.fpkm.xls', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GeneExp/{sample}/{sample}.Bowtie2Gene.MapReadsStat.xls', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GeneExp/{sample}/{sample}.ReadsRandom.png', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/GeneExp/{sample}/{sample}.SeqSaturation.png', outdir=outdir, sample=SAMPLES),
        # expand('{outdir}/GO_Hypergeometric/tmp_file/{repcompare}.DEseq2_Method.glist',outdir=outdir,repcompare=repcompare),
        # expand('{outdir}/GO_Hypergeometric/tmp_file/{repcompare}.DEseq2_Method.glist2',outdir=outdir,repcompare=repcompare),
        # expand('{outdir}/GO_Hypergeometric/GO/{repcompare}.DEseq2_Method_BP_EnrichResult.sorted.padjust.xls',outdir=outdir,repcompare=repcompare),
        # expand('{outdir}/GO_Hypergeometric/GO/{repcompare}.DEseq2_Method_CC_EnrichResult.sorted.padjust.xls',outdir=outdir,repcompare=repcompare),
        # expand('{outdir}/GO_Hypergeometric/GO/{repcompare}.DEseq2_Method_MF_EnrichResult.sorted.padjust.xls',outdir=outdir,repcompare=repcompare)

#可以通过flag来控制是否需要运行某一个模块
if filter_step == 1:
    include: ruledir+'/Filter.snakemake.py'
if alignment_step == 1:
    include: ruledir+'/Alignment.snakemake.py'
if geneexp_step == 1:
    include: ruledir+'/GeneExp.snakemake.py'
if genediff_step == 1 and rep == 1 :
    include: ruledir+'/GeneDiff.snakemake.py'
if genediff_step == 1 and rep == 0 :
    include: ruledir + '/GeneDiff.norepeat.snakemake.py'
if go_step == 1 :
    include: ruledir + '/GO.snakemake.py'