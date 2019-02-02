#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: rule.test.snakemake.py.py
@time: 2019/01/04
'''

rule 01:
    input:
        'F1.sh'
    output:
        'F1.sh.sign'
    run：
        shell('sh F1.sh')
rule 02:
    input:
        'F1.sh.sign'
    output:
        'F2.sh.sign'
    run:
        shell('sh F2.sh')

#如果是已经处理过的sh，那么只需要根据这个rule写后面的依赖关系即可
rule 03:
    input:
        'F1.sh.sign'
    output:
        'F3.sh.sign'
    run:
        shell("sh F3.sh")
