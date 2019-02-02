#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Target.snakemake.py
@time: 2019/01/11
'''


outdir='/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3'
SAMPLES=config["samples"]


rule all:
    input:
        expand('{outdir}/Target_miRanda/{sample}.circBase-miRBase.xls',outdir=outdir,sample=SAMPLES)

#only for hsa and mus
rule target:
    input:
        '/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/database/miRBase-circBase/hsa/hsa.miRanda.xls'
    output:
        expand('{outdir}/Target_miRanda/{sample}.circBase-miRBase.xls',outdir=outdir,sample=SAMPLES)
    params:
        sb='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Target/',
        od=expand('{outdir}',outdir=outdir)
    shell:
        'perl {params.sb}/extract_circRNA-miRBase.pl {input} {params.od}/Annotation {params.od}/Target_miRanda'