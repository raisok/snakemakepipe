#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GeneExp.snakemake_md.py.py
@time: 2018/12/12
'''


# configfile: "rnaseq.yaml"

SAMPLES = config["samples"]
soft=config["software"]
outdir=config["outdir"]
sb=config["software"]["scriptbin"]
gene=config["database"]["gene"][0]
gene2tr=config["database"]["gene"][1]

# rule all:
#     input:
#         expand('{outdir}/GeneExp/{sample}/{sample}.Map2GeneStat.xls',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.isoforms.fpkm.xls',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.Bowtie2Gene.MapReadsStat.xls',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.ReadsRandom.png',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.SeqSaturation.png',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/GeneExp/{sample}/{sample}.ReadsCoverage.png',outdir=outdir,sample=SAMPLES),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/AllSamples.GeneExpression.FPKM.xls',outdir=outdir),
#         #expand('{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/GeneExpressionSummary.xls',outdir=outdir),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/CorrelationHeatmap/AllSamples.CorrelationHeatmap.png',outdir=outdir),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/HclusterTree/AllSamples.HclustTree.png',outdir=outdir),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Density.png',outdir=outdir),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Expression_distribution.png',outdir=outdir),
#         expand('{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Boxplot.png',outdir=outdir)

rule buildindex:
    input:
        gene2tr=gene2tr,
        cds=gene
    output:
        '{outdir}/GeneExp/rsem-build/refMrna.fa',
        '{outdir}/GeneExp/rsem-build/gene2tr.txt',
        '{outdir}/GeneExp/TranscriptLength.bed',
        '{outdir}/GeneExp/TranscriptLength.txt'
    params:
        bowtie2 = soft["bowtie2"][0],
        rsem = soft["rsem"][0],
        od = outdir + "/GeneExp",
        scriptbin = sb[0],
        libpath="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/RNA_lib/"
    log:
        e = outdir + "/log/buildindex.e",
        o = outdir + "/log/buildindex.o",
        time = outdir + "/time/buildindex.time"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("cat {input.cds} >{output[0]}; \
               cat {input.gene2tr} >{output[1]}; \
              export LD_LIBRARY_PATH={params.libpath}:$LD_LIBRARY_PATH; \
              {params.rsem}/rsem-prepare-reference {output[0]} {output[0]} \
              --bowtie2 --bowtie2-path {params.bowtie2} --transcript-to-gene-map {output[1]}; \
              perl {params.scriptbin}/GeneExp/fastaDeal.pl -attr id:len {output[0]} >{output[3]}; \
              awk '{{print $1\"\\t1\\t\"$2}}' {output[3]} >{output[2]}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule geneexp:
    input:
        fq1 ='{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz',
        fq2 ='{outdir}/Filter_SOAPnuke/{sample}/{sample}.clean.1.fq.gz',
        refmrna='{outdir}/GeneExp/rsem-build/refMrna.fa',
        tranlen='{outdir}/GeneExp/TranscriptLength.txt',
        tranbed='{outdir}/GeneExp/TranscriptLength.bed'
    output:
        '{outdir}/GeneExp/{sample}/{sample}.bam',
        '{outdir}/GeneExp/{sample}/{sample}.Map2GeneStat.xls',
        '{outdir}/GeneExp/{sample}/{sample}.genes.results',
        '{outdir}/GeneExp/{sample}/{sample}.isoforms.results',
        '{outdir}/GeneExp/{sample}/{sample}.gene.fpkm.xls',
        '{outdir}/GeneExp/{sample}/{sample}.isoforms.fpkm.xls',
        '{outdir}/GeneExp/{sample}/{sample}.transcript.bam',
        '{outdir}/GeneExp/{sample}/{sample}.transcript.sorted.bam'
    params:
        bowtie2=soft["bowtie2"][0],
        samtools=soft["samtools"][0],
        rsem = soft["rsem"][0],
        para=soft["bowtie2"][1],
        od=outdir+"/GeneExp",
        scriptbin=sb[0],
        libpath="/ldfssz1/ST_BIGDATA/PMO/SOFTWARE/RNA_SoftWare/RNA_lib/"
    log:
        e=outdir+"/log/geneexp_{sample}.e",
        o=outdir+"/log/geneexp_{sample}.o",
        time=outdir+"/time/geneexp_{sample}.time"
    threads:24
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell(" export LD_LIBRARY_PATH={params.libpath}:$LD_LIBRARY_PATH; \
              {params.bowtie2}/bowtie2 {params.para} -1 {input.fq1} -2 {input.fq2} -x {input.refmrna} 2>{output[1]}| \
            {params.samtools} view -b -S -o {output[0]} - ; \
            {params.rsem}/rsem-calculate-expression --paired-end -p 8 --bam {output[0]} {input.refmrna} {params.od}/{wildcards.sample}/{wildcards.sample}; \
            awk '{{if($7!=0.00)print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' {output[2]} |grep -v '^ERCC' >{output[4]}; \
            awk '{{if($7!=0.00)print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"$7}}' {output[3]} |grep -v '^ERCC' >{output[5]}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule geneexp_stat:
    input:
        bam ='{outdir}/GeneExp/{sample}/{sample}.bam',
        gene2tr=gene2tr
    output:
        '{outdir}/GeneExp/{sample}/{sample}.Bowtie2Gene.MapReadsStat.xls'
    log:
        e=outdir+"/log/gene_stat.e",
        o=outdir+"/log/gene_stat.o",
        time=outdir+"/time/gene_stat.time"
    params:
        scriptbin=sb[0],
        samtools=soft["samtools"][0],
        od=outdir+"/GeneExp",
        PEorSE="PE"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.scriptbin}/GeneExp/BowtieMapStat.pl -bam {input.bam} -key \
        {params.od}/{wildcards.sample}/{wildcards.sample}.Bowtie2Gene --seqType {params.PEorSE} -samtools {params.samtools} \
         -gene2tr {input.gene2tr}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule geneexp_random:
    input:
        bam = '{outdir}/GeneExp/{sample}/{sample}.bam',
        len = '{outdir}/GeneExp/TranscriptLength.txt',
        gene2tr = gene2tr
    output:
        '{outdir}/GeneExp/{sample}/{sample}.ReadsRandom.pdf',
        '{outdir}/GeneExp/{sample}/{sample}.ReadsRandom.png'
    log:
        e = outdir + "/log/gene_random.e",
        o = outdir + "/log/gene_random.o",
        time = outdir + "/time/gene_random.time"
    params:
        scriptbin = sb[0],
        samtools = soft["samtools"][0],
        od=outdir+"/GeneExp",
        PEorSE = "PE"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.scriptbin}/GeneExp/ReadsRandomInGene_forBowtie_ggplot2.pl -len {input.len} -bam {input.bam} \
        -seqType {params.PEorSE} -prefix {params.od}/{wildcards.sample}/{wildcards.sample}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule geneexp_cov:
    input:
        bam = '{outdir}/GeneExp/{sample}/{sample}.transcript.sorted.bam',
        len = '{outdir}/GeneExp/TranscriptLength.txt',
        bed = '{outdir}/GeneExp/TranscriptLength.bed'
    output:
        '{outdir}/GeneExp/{sample}/{sample}.depth.txt',
        '{outdir}/GeneExp/{sample}/{sample}.ReadsCoverage.pdf',
        '{outdir}/GeneExp/{sample}/{sample}.ReadsCoverage.png'
    log:
        e = outdir + "/log/gene_cov.e",
        o = outdir + "/log/gene_cov.o",
        time = outdir + "/time/gene_cov.time"
    params:
        scriptbin = sb[0],
        samtools = soft["samtools"][0],
        od=outdir+"/GeneExp",
        PEorSE = "PE"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("{params.samtools} depth -b {input.bed} {input.bam} > {output[0]}; \
              perl {params.scriptbin}/GeneExp/cover_stat.pl {output[0]} {input.len} {wildcards.sample} {params.od}/{wildcards.sample}")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

rule geneexp_seq:
    input:
        bam = '{outdir}/GeneExp/{sample}/{sample}.transcript.bam',
        gene2tr = gene2tr
    output:
        '{outdir}/GeneExp/{sample}/{sample}.SeqSaturation.png'
    log:
        e = outdir + "/log/gene_seq.e",
        o = outdir + "/log/gene_seq.o",
        time = outdir + "/time/gene_seq.time"
    params:
        scriptbin = sb[0],
        samtools = soft["samtools"][0],
        od=outdir+"/GeneExp",
        PEorSE = "PE"
    run:
        shell("echo start [ `date +%F` `date +%T` ] >{log.time}")
        shell("perl {params.scriptbin}/GeneExp/SeqSaturation_Bowtie.pl -bam {input.bam} -seqType {params.PEorSE} -gene2tr \
                {input.gene2tr} -prefix {params.od}/{wildcards.sample}/{wildcards.sample} -cutoff 1 ")
        shell("echo finish [ `date +%F` `date +%T` ] >>{log.time}")

#暂时未添加PCA和Venn
rule gene_statistic:
    input:
        indir='{outdir}/GeneExp/'
    output:
        '{outdir}/GeneExp/AllSamples.GeneExpression.FPKM.xls',
        '{outdir}/GeneExp/AllSamples.TranscriptExpression.FPKM.xls',
        '{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/AllSamples.GeneExpression.FPKM.xls',
        '{outdir}/BGI_result/Quantify/GeneExpression/CorrelationHeatmap/AllSamples.CorrelationHeatmap.png',
        '{outdir}/BGI_result/Quantify/GeneExpression/HclusterTree/AllSamples.HclustTree.png',
        '{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Density.png',
        '{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Expression_distribution.png',
        '{outdir}/BGI_result/Quantify/GeneExpression/GeneExpression/Boxplot.png'
    params:
        scriptbin = sb[0],
        samtools = soft["samtools"][0],
        od = outdir + "/BGI_result/Quantify/GeneExpression/",
        PEorSE = "PE"
    run:
        odtmp=os.path.dirname(output[2])
        shell("perl {params.scriptbin}/GeneExp/AllGeneStat.pl {input.indir} {output[0]}; \
              perl {params.scriptbin}/GeneExp/AllTranscriptStat.pl {input.indir} {output[1]}; \
              cp {output[0]} {output[2]}; \
              perl {params.scriptbin}/GeneExp/drawCorrelationHeatmap.pl -indir {odtmp} -outdir {params.od}/CorrelationHeatmap; \
              perl {params.scriptbin}/GeneExp/drawHclustTree.pl -indir {odtmp} -outdir {params.od}/HclusterTree; \
              perl {params.scriptbin}/GeneExp/drawbox_density-plot.pl -i {odtmp} -o {odtmp}; \
              perl {params.scriptbin}/GeneExp/drawstacked_bar.pl -indir {odtmp} -outdir {odtmp}")