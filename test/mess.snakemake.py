#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: mess.snakemake_md.py.py
@time: 2019/01/02
'''

os.system('touch a.input')

rule NAME:
    input:
        "a.input"
    output:
        "b.output"
    threads: 1
    message: "executing cat command with {threads} threads on {input} file."
    shell: "cat {input} > {output}"