#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Rm_rRNA.Filter.snakemake.py
@time: 2019/01/03
'''
import os

SAMPLES = config["samples"]
soft=config["software"]
sb=config["software"]["scriptbin"]

def getfq1(samples):
    for i in samples:
        return samples[i][0]
def getfq2(samples):
    for i in samples:
        return samples[i][1]

rRNAfa="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Human_rRNA_NCBI.fa"
outdir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3"

# rule all:
#     input:
#         expand('{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved_1.fq.gz',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved_2.fq.gz',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.base.png', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.qual.png', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.RawReadsClass.png', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.filter.stat.xls', outdir=outdir, sample=SAMPLES),
#         expand('{outdir}/Filter_SOAPnuke/FilterSummary.xls', outdir=outdir)

rule builder:
    input:
        rRNAfa
    output:
        '{outdir}/Rm_rRNA/index/Human_rRNA_NCBI.fa.index.rev.pac'
    params:
        soap="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Rm_rRNA/bin/Soap/2bwt-builder",
    run:
        od=os.path.dirname(output[0])
        shell("ln -s {input} {od}/Human_rRNA_NCBI.fa; \
        {params.soap} {od}/Human_rRNA_NCBI.fa")

rule Rm_rRNA:
    input:
        fq1=lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2=lambda wildcards: SAMPLES[wildcards.sample][1],
        index=rules.builder.output
    output:
        spfq1='{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNA.PESoap.gz',
        spfq2='{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNA.PESoapSingle.gz',
        rmrnafq1='{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved_1.fq.gz',
        rmrnafq2='{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved_2.fq.gz'
    params:
        soap="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Rm_rRNA/bin/Soap/soap_mm_gz",
        para="-m 0 -x 1000 -s 28 -l 32 -v 5 -r 2 -p 3",
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Rm_rRNA/bin",
        sb2=sb[0]+"/Filter/",
        prefix='{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved',
        od='{outdir}/Rm_rRNA/{sample}/{sample}_Lib',
        fqcheck='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Rm_rRNA/bin/fqcheck'
    run:
        idx=os.path.dirname(input.index[0])+"/"+os.path.basename(input.index[0]).split('.')[0]+".fa.index"
        shell("{params.soap} -a {input.fq1} -b {input.fq2} -D {idx} {params.para} -o {output.spfq1} -2 {output.spfq2}; \
              perl {params.sb}/rRNAFilter.pl -fq {input.fq1},{input.fq2} -soap {output.spfq1},{output.spfq2} -output {params.prefix}; \
              {params.fqcheck} -r {params.od}/{wildcards.sample}_Lib_rRNAremoved_1.fq.gz -c {params.od}/1.fqcheck; \
              {params.fqcheck} -r {params.od}/{wildcards.sample}_Lib_rRNAremoved_2.fq.gz -c {params.od}/2.fqcheck; \
              perl {params.sb2}/fqcheck_distribute.pl {params.od}/1.fqcheck {params.od}/2.fqcheck -o {params.prefix}. ;")

rule rRNA_stat:
    input:
        '{outdir}/Rm_rRNA/{sample}/{sample}_Lib/{sample}_Lib_rRNAremoved.rrna.stat'
    output:
        '{outdir}/Rm_rRNA/{sample}/{sample}_rRNAremoved.rrna.stat'
    params:
        sb="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Rm_rRNA/bin"
    run:
        shell("perl  {params.sb}/stat.pl -in {input} -out {output}")

rule filter_fastq:
    input:
        fq1=rules.Rm_rRNA.output.rmrnafq1,
        fq2=rules.Rm_rRNA.output.rmrnafq2
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
              perl {params.scriptbin}/soapnuke_stat.pl {params.od}/{wildcards.sample}/  >{wildcards.sample}.filter.stat.xls; \
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
