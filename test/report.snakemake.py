#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: report.snakemake_md.py
@time: 2018/12/26
'''

rule report:
    input:
        "rnaseq.yaml"
    output:
        "report.html"
    run:
        from snakemake_md.utils import report
        with open(input[0]) as test:
            n_calls = sum (1 for l in test if not l.startswith('#'))
        report(
            '''
            An example variant calling workflow
            ===================================
            
            Reads were mapped to the Yeast reference genome ang variants were called jointly with SAMtools/BCFtools.
            
            The resulted in {n_calls} variants (see Table T1_).
            ''',output[0],T1=input[0]
        )