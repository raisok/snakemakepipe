#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: KEGG.enrich.snakemake.py
@time: 2019/01/11
'''

outdir='/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example/result3'

method=config['diff_para']['methods']
repcompare=config['diff_para']['repcompare']
norepcompare=config['diff_para']['norepcompare']
group=config['diff_para']['group']

rule all:
    input:
        expand('{outdir}/Pathway_Hypergeometric/Pathway',outdir=outdir)


rule keggenrichment:
    input:
        expand('{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExpFilter.xls',outdir=outdir,repcompare=repcompare),
        expand('{outdir}/GeneDiffExp_Allin/DEGseq/{norepcompare}.DEGseq_Method.DiffExpFilter.xls',outdir=outdir,norepcompare=norepcompare)
    output:
        '{outdir}/Pathway_Hypergeometric/Pathway'
    params:
        keggFa='/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/animal.fa',
        koMap='/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/komap/animal_ko_map.tab',
        mapTitle='/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map_title.tab',
        mapDir='/ifs4/BC_PUB/biosoft/db/Pub/kegg/RNA/84.0/map',
        outdir='{outdir}/Pathway_Hypergeometric/Pathway',
        tmpdir='{outdir}/Pathway_Hypergeometric/tmp_file',
        bg='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19.ko',
        sb='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Enrichment/',
        tmp='{outdir}/Pathway_Hypergeometric/tmp'
    shell:
        "mkdir -p {params.tmp};mkdir -p {params.tmpdir};cp {input} {params.tmp}/; \
keggFa={params.keggFa};koMap={params.koMap};mapTitle={params.mapTitle};mapDir={params.mapDir};outdir={params.outdir};tmpdir={params.tmpdir};bg={params.bg}; \
for i in `ls {params.tmp}/*.DiffExpFilter.xls`;do keyname=`echo $i|sed 's/.*\/\(.*\).DiffExpFilter.xls/\\1/'`; \
awk '{{print $1\" \"$5}}' $i >$tmpdir/$keyname.glist; \
perl {params.sb}/getKO.pl -glist $tmpdir/$keyname.glist -bg $bg -outdir $tmpdir; \
perl {params.sb}/pathfind.pl -kegg $keggFa -komap $koMap -maptitle $mapTitle -fg $tmpdir/$keyname.ko -bg $bg -output $outdir/$keyname.path; \
mkdir -p $outdir/$keyname\_map; \
perl {params.sb}/keggMap.pl -ko $tmpdir/$keyname.ko -diff $tmpdir/$keyname.glist -komap $koMap -mapdir $mapDir -outdir $outdir/$keyname\_map; \
awk '{{if($5>0)print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\"\\tdown\"}}' $i >$tmpdir/$keyname.glist.temp; \
perl {params.sb}/drawKEGG.pl -path $outdir/$keyname.path -outprefix $outdir/$keyname -idCol 6 -level1Col 7 -level2Col 8 -geneCol 9 -list $tmpdir/$keyname.glist.temp; \
perl {params.sb}/Pathway_gene2symbol.pl -symbol /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/human.gene2symbol.txt -G2K $outdir/$keyname.Gene2KEGG.xls -K2G $outdir/$keyname.KEGG2Gene.xls -outprefix $outdir/$keyname; \
done; \
perl {params.sb}/genPathHTML.pl -indir $outdir; \
perl {params.sb}/pathway_enrichFigure.pl $outdir \
"