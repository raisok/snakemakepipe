#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: cpfq.snakemake_md.py
@time: 2018/12/28
'''

SAMPLES = config["samples"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]

rule all:
    input:
        expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz', outdir=outdir,sample=SAMPLES),
        expand('{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz', outdir=outdir,sample=SAMPLES),

rule cp_fq:
    input:
        fq1 = lambda wildcards: SAMPLES[wildcards.sample][0],
        fq2 = lambda wildcards: SAMPLES[wildcards.sample][1]
    output:
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz',
        '{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.2.fq.gz'
    run:
        shell('ln -s {input.fq1} {output[0]};ln -s {input.fq2} {output[1]};')