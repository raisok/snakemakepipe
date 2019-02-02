#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: test.snakemake_md.py
@time: 2018/12/25
'''


rule all:
    input:
        "a.flag2",
        "b.flag",
        "c.flag"

rule test:
    output:
        touch('a.flag2')
    script:
        "test.py"


rule test2:
    output:
        touch('b.flag')
    shell:
        "echo this is test2"


rule test3:
    output:
        touch('c.flag')
    log:
        "test3.log"
    run:
        shell("perl test.pl ")
        shell("perl -e 'print \"this is test3\"; ' >>{log}")