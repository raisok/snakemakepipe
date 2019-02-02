#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import time
import os

'''
@author: yueyao
@file: Alignment.snakemake_md.py.py
@time: 2018/12/07
'''

# configfile: "rnaseq.yaml"

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]
gn=config["database"]["genome"][1]
gnidx=config["database"]["genome"][0]

# rule all:
#     input:
#         expand('{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam.bai',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GenomeMapping_HISAT/GenomeMappingSummary.xls', outdir=outdir)

rule alignment:
    input:
        fq1 ='{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz',
        fq2 ='{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz'
    output:
        '{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam',
        '{outdir}/GenomeMapping_HISAT/{sample}/{sample}.AddRG.Reorder.Sort.bam.bai',
        '{outdir}/GenomeMapping_HISAT/{sample}/{sample}.Map2GenomeStat.xls'
    params:
        hisat2=soft["hisat2"][0],
        java=soft["java"][0],
        samtools=soft["samtools"][0],
        para=soft["hisat2"][1],
        picard=soft["picard"][0],
        od=outdir+"/GenomeMapping_HISAT",
        genome=gnidx,
        scriptbin=sb[0]+'GenomeMapping_HISAT/',
        gid=gn
    log:
        e=outdir+"/log/alignment_{sample}.e",
        o=outdir+"/log/alignment_{sample}.o",
        time=outdir+"/time/alignment_{sample}.time"
    threads:24
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("cd {params.od}/{wildcards.sample}; \
        {params.hisat2} {params.para} -p {threads} -1 {input.fq1} -2 {input.fq2} -x {params.genome} 2>{wildcards.sample}.Map2GenomeStat.xls| \
        {params.samtools} view -b -S -o {wildcards.sample}.bam -; \
        if [ ! -d java_tmp ];then mkdir -p java_tmp;fi; \
              {params.java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {params.picard}/AddOrReplaceReadGroups.jar I={wildcards.sample}.bam O={wildcards.sample}.AddRG.bam RGID={wildcards.sample} RGLB={wildcards.sample}_library RGPL=illumina RGPU=machine RGSM={wildcards.sample} VALIDATION_STRINGENCY=SILENT; \
              {params.java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {params.picard}/ReorderSam.jar I={wildcards.sample}.AddRG.bam O={wildcards.sample}.AddRG.Reorder.bam R={params.gid} VALIDATION_STRINGENCY=SILENT; \
              {params.java} -Xmx4G -Djava.io.tmpdir=java_tmp -jar {params.picard}/SortSam.jar I={wildcards.sample}.AddRG.Reorder.bam O={wildcards.sample}.AddRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT; \
              {params.samtools} index {output[0]} \
              1>>{log.o} 2>>{log.e}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule alignment_stat:
    input:
        expand('{outdir}/GenomeMapping_HISAT/{sample}/{sample}.Map2GenomeStat.xls', outdir=outdir,sample=SAMPLES)
    output:
        outdir+"/GenomeMapping_HISAT/GenomeMappingSummary.xls"
    log:
        e=outdir+"/log/alignment_stat.e",
        o=outdir+"/log/alignment_stat.o",
        time=outdir+"/time/alignment_stat.time"
    params:
        scriptbin=sb[0]+'/Alignment/',
        peorse="PE"
    run:
        tmpindir = os.path.dirname(input[0])
        indir = os.path.dirname(tmpindir)
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.scriptbin}/MapStat.pl -indir  {indir} -output {output} -type {params.peorse}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")