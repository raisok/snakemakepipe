#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Denovo.stLFR.snakemake.py.py
@time: 2018/12/12
'''

SAMPLES=config["samples"]
soft=config["software"]
outdir=config["outdir"]



rule all:
    input:
        expand('{outdir}/{sample}/Denovo_Trinity/Unigene.fa', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/{sample}/Denovo_Trinity/All-Unigene.fa', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/{sample}/{sample}.Map2TrinityStat.xls', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/{sample}/{sample}.bam', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/Map2Trinity.xls', outdir=outdir),
        expand('{outdir}/assembly_trinity.stat.xls',outdir=outdir),
        expand('{outdir}/assembly_unigene.stat.xls',outdir=outdir)

rule trinity_denovo:
    input:
        leftfq=lambda wildcards: SAMPLES[wildcards.sample][0],
        rightfq=lambda wildcards: SAMPLES[wildcards.sample][1]
    output:
        '{outdir}/{sample}/Denovo_Trinity/Unigene.fa',
        '{outdir}/{sample}/Denovo_Trinity/All-Unigene.fa'
    params:
        od=outdir,
        trinity=soft["trinity"][0],
        para=soft["trinity"][1],
        python='/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/bin/python',
        sb='/ldfssz1/ST_BIGDATA/USER/yueyao/bin/'
    log:
        e = outdir + "/log/trinity_{sample}.e",
        o = outdir + "/log/trinity_{sample}.o",
        time = outdir + "/time/trinity_{sample}.time"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("export PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/jre1.8.0_45/bin/:/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/bowtie2-2.2.5/:/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/samtools-0.1.19/$PATH; \
            export LD_LIBRARY_PATH=/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/lib/:$LD_LIBRARY_PATH; \
            mkdir -p {params.od}/{wildcards.sample}/Denovo_Trinity/; \
            cd /{params.od}/{wildcards.sample}/Denovo_Trinity/; \
            {params.trinity} {params.para}  \
            --left {input.leftfq} --right {input.rightfq} --output {params.od}/{wildcards.sample}/Denovo_Trinity/Trinity; \
            ln -s {params.od}/{wildcards.sample}/Denovo_Trinity/Trinity/Trinity.fasta Unigene.fa; \
            {params.python} \
            {params.sb}/extract_longest_isform.py Unigene.fa All ; \
            perl {params.sb}/fishInWinter.pl -bf table -ff fasta All_Unigene_id.txt Unigene.fa >Temp.Unigene.fa; \
            perl -lane 'if(/(>.*c\d+_g\d+)(_i\d+)\slen/){{print $1}}else{{print}}' Temp.Unigene.fa >All-Unigene.fa; \
            perl {params.sb}/get_Trinity_gene_to_trans_map.pl Unigene.fa >All-Unigene.gene2mark; \
            rm Temp.Unigene.fa")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule bowtie2unigene:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2 = lambda wildcards: SAMPLES[wildcards.sample][1],
        unigene='{outdir}/{sample}/Denovo_Trinity/Unigene.fa',
    output:
        '{outdir}/{sample}/{sample}.Map2TrinityStat.xls',
        '{outdir}/{sample}/{sample}.bam'
    params:
        od=outdir,
        bowtie2=soft["bowtie2"][0]+"bowtie2",
        para=soft["bowtie2"][1],
        samtools=soft["samtools"]
    log:
        time = outdir + "/time/bowtie2_{sample}.time"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("mkdir -p {params.od}/{wildcards.sample}/index; \
              {params.bowtie2}/bowtie2-build \
              {input.unigene} \
              {params.od}/{wildcards.sample}/index/refMrna.fa; \
              {params.bowtie2}/bowtie2 \
              {params.para} -p {threads}\
              -x {params.od}/{wildcards.sample}/index/refMrna.fa \
              -1 {input.fq1} \
              -2 {input.fq2} \
              2>{output[0]} | \
              {params.samtools} view -S -b -o {output[1]} - ")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule bowtiestat:
    input:
        expand('{outdir}/{sample}/{sample}.Map2TrinityStat.xls',outdir=outdir,sample=SAMPLES)
    output:
        '{outdir}/Map2Trinity.xls'
    params:
        outdir=outdir,
        sb='/ldfssz1/ST_BIGDATA/USER/yueyao/bin/',
        peorse='pe'
    log:
        time = outdir + "/time/bowtiestat_{wildcards.sample}.time"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.sb}/MapStat.pl -indir {params.outdir} -output {output[0]} -tpye {params.peorse}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule trinitystat:
    input:
        expand('{outdir}/{sample}/Denovo_Trinity/Unigene.fa', outdir=outdir, sample=SAMPLES),
        expand('{outdir}/{sample}/Denovo_Trinity/All-Unigene.fa', outdir=outdir, sample=SAMPLES),
    output:
        '{outdir}/assembly_trinity.stat.xls',
        '{outdir}/assembly_unigene.stat.xls'
    params:
        id=outdir,
        sb='/ldfssz1/ST_BIGDATA/USER/yueyao/bin/'
    log:
        time = outdir + "/time/trinitystat_{wildcards.sample}.time"
    threads: 1
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.sb}/assembly_stat2.pl -indir {params.id} -outdir {params.id}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")