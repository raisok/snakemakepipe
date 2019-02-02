#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: Prediction.snakemake.py
@time: 2019/01/09
'''

include: 'Rm_rRNA.Filter.snakemake.py'

SAMPLES = config["samples"]
soft=config["software"]
sb=config["software"]["scriptbin"]
outdir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3"

# rule all:
#     input:
#         expand('{outdir}/Prediction/{sample}/{sample}_ciri_output.xls',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/Prediction/{sample}/{sample}_findcirc_output.xls',outdir=outdir,sample=SAMPLES)

# CIRI prediction v2.0.5
rule ciri_predict:
    input:
        fq1=rules.filter_fastq.output[0],
        fq2=rules.filter_fastq.output[1]
    output:
        '{outdir}/Prediction/{sample}/{sample}_ciri_output.xls'
    params:
        bwa="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/envs/rnaseq/bin/bwa",
        genomeid="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GenomeBwaIndex/chrALL.fa",
        para="mem -T 19",
        perl="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/envs/rnaseq/bin/perl",
        ciri="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Prediction/../software/CIRI_v2.0.5/CIRI_v2.0.5.pl"
    run:
        od=os.path.dirname(output[0])
        shell('cd {od};{params.bwa} {params.para} {params.genomeid} {input.fq1} {input.fq2} 1>{wildcards.sample}_bwa.sam 2> {wildcards.sample}_bwa.log; \
              {params.perl} {params.ciri} -I {wildcards.sample}_bwa.sam -O {output} -F {params.genomeid} ')

# find_circ prediction
rule findcirc_prediction:
    input:
        fq1 = rules.filter_fastq.output[0],
        fq2 = rules.filter_fastq.output[1]
    output:
        '{outdir}/Prediction/{sample}/{sample}_findcirc_output.xls',
        '{outdir}/Prediction/{sample}/{sample}_bowtie2.log'
    params:
        bowtie2="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/envs/rnaseq/bin/bowtie2",
        para="--very-sensitive --score-min=C,-15,0 --mm",
        genomeid="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/prepare/GenomeBowtie2Index/chrALL.fa",
        path="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNAref/software/anaconda3/envs/rnaseq/bin",
        findcirc="/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Prediction/../software/find_circ/trunk/"
    run:
        od=os.path.dirname(output[0])
        shell('cd {od}; \
              export PATH={params.path}:$PATH; \
              {params.bowtie2} {params.para} -x {params.genomeid} -1 {input.fq1} -2 {input.fq2} 2>{output[1]} | samtools view -hbuS - | samtools sort - {wildcards.sample}_bowtie2.sort; \
               samtools view -hf 4 {wildcards.sample}_bowtie2.sort.bam | samtools view -Sb - > unmapped_{wildcards.sample}.bam; \
              {params.findcirc}/unmapped2anchors.py unmapped_{wildcards.sample}.bam | gzip > {wildcards.sample}_anchors.fq.gz; \
              bowtie2 --score-min=C,-15,0 --reorder --mm -x {params.genomeid} -U {wildcards.sample}_anchors.fq.gz | {params.findcirc}/find_circ.py --genome={params.genomeid} --name={wildcards.sample} --stats=./stat.txt --reads=./spliced_reads.fa > ./splice_sites.bed; \
              grep \'^#\' ./splice_sites.bed | uniq > {output[0]}; \
              grep CIRCULAR ./splice_sites.bed | grep -v chrM | awk \'$5>=2\' | grep UNAMBIGUOUS_BP | grep ANCHOR_UNIQUE | {params.findcirc}/maxlength.py 100000 >> {output[0]}')