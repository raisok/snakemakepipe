#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import time
import os

'''
@author: yueyao
@file: Filter.snakemake_md.py
@time: 2018/12/05
'''

# configfile:"rnaseq.yaml"
# configfile:"test.json"
#
# snakemake_md -config "rnaseq.yaml"
# snakemake_md -config "test.json"

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]

# rule all:
#     input:
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.base.png', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.qual.png', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.RawReadsClass.png', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.filter.stat.xls', outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/FilterSummary.xls', outdir=outdir)

rule filter_fastq:
    input:
        fq1=lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2 =lambda wildcards: SAMPLES[wildcards.sample][1]
    output:
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.base.png',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.qual.png',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.RawReadsClass.png',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.filter.stat.xls'
    params:
        soapnuke=soft["soapnuke"][0],
        fqcheck=soft["fqcheck"][0],
        para=soft["soapnuke"][1],
        od=outdir+"/Filter_SOAPnuke",
        scriptbin=sb[0]+"/Filter/"
    log:
        e=outdir+"/log/filter_{sample}.e",
        o=outdir+"/log/filter_{sample}.o",
        time=outdir+"/time/filter_{sample}.time"
    # threads:1
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("{params.soapnuke} filter {params.para} -1 {input.fq1} -2 {input.fq2} -o {params.od}/{wildcards.sample} -C {wildcards.sample}.clean.1.fq.gz -D {wildcards.sample}.clean.2.fq.gz -R {wildcards.sample}.raw.1.fq.gz -W {wildcards.sample}.raw.2.fq.gz 1>{log.o} 2>{log.e}")
        shell("cd {params.od}/{wildcards.sample}; \
              {params.fqcheck} -r {wildcards.sample}.clean.1.fq.gz -c {wildcards.sample}_1.fqcheck; \
              {params.fqcheck} -r {wildcards.sample}.clean.2.fq.gz -c {wildcards.sample}_2.fqcheck; \
              {params.fqcheck} -r {wildcards.sample}.raw.1.fq.gz -c {wildcards.sample}_1.rawdata.fqcheck; \
              {params.fqcheck} -r {wildcards.sample}.raw.2.fq.gz -c {wildcards.sample}_1.rawdata.fqcheck; \
              perl {params.scriptbin}/fqcheck_distribute.pl {wildcards.sample}_1.fqcheck {wildcards.sample}_2.fqcheck -o {wildcards.sample}. ; \
              perl {params.scriptbin}/soapnuke_stat.pl {params.od}/{wildcards.sample} >{wildcards.sample}.filter.stat.xls; \
              perl {params.scriptbin}/drawPizza.pl -infile {params.od}/{wildcards.sample}/{wildcards.sample}.filter.stat.xls -outdir {params.od}/{wildcards.sample} \
              1>>{log.o} 2>>{log.e}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule filter_stat:
    input:
        expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.filter.stat.xls', outdir=outdir,sample=SAMPLES)
    output:
        outdir+"/Filter_SOAPnuke/FilterSummary.xls"
    log:
        e=outdir+"/log/filter_stat.e",
        o=outdir+"/log/filter_stat.o",
        time=outdir+"/time/filter_stat.time"
    params:
        scriptbin=sb[0]+"/Filter/"
    run:
        tmpindir=os.path.dirname(input[0])
        indir=os.path.dirname(tmpindir)
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.scriptbin}/filter_stat.pl -indir {indir} -output {output}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")