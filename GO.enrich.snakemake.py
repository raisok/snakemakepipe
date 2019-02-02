#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: GO.enrich.snakemake.py
@time: 2019/01/11
'''

rule goenrichment:
    input:
        '{outdir}/GeneDiffExp_Allin/DEseq2/{repcompare}.DEseq2_Method.DiffExpFilter.xls'
    output:
        ''
    params:
        goclass='/ifs4/BC_PUB/biosoft/db/Pub/go/RNA/20171220/go.class',
        outdir='{outdir}/GO_Hypergeometric/GO',
        tmpdir='{outdir}/GO_Hypergeometric/tmp_file',
        prefix='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg79/hg19',
        gene2tr='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/Annotation_kegg76/refMrna.fa.gene2mark',
        path='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAdenovo/RNA_RNAdenovo_2015a/software/perl-V5/bin',
        sb='/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_circRNA/RNA_circRNA_2017a/Enrichment/'
    shell:
        "export PATH={params.path}:$PATH; \
         perl {params.sb}/merge_gene2tr.pl {params.gene2tr} {params.tmpdir}/gene2tr; \
         gene2tr={params.tmpdir}/gene2tr; \
         perl {params.sb}/merge_go.pl $gene2tr {params.prefix} {params.tmpdir}/species; \
         prefix={params.tmpdir}/species; \
         outdir={params.outdir}; \
         tmpdir={params.tmpdir}; \
         for i in `ls {}/*.GeneDiff*Filter.xls`;do keyname=`echo $i|sed 's/.*\/\(.*\).GeneDiff.*Filter.xls/\\1/`; \
         'awk '{{if($5>0)print $1\"\\t\"$5\"\\tup\";else print $1\"\\t\"$5\" down"}}' $i >$tmpdir/$keyname.glist; \
         perl {params.sb}/drawGO.pl -list $tmpdir/$keyname.glist -goclass $goclass -goprefix $prefix -outprefix $outdir/$keyname; done; \
         perl  go.pl -gldir $tmpdir -sdir `dirname $prefix` -species `basename $prefix` -outdir $outdir; \
         perl topGO.pl -gldir $tmpdir -godir $outdir -prefix $prefix -list $gene2tr -outdir $outdir \
         "
