#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''
@author: yueyao
@file: circos.snakemake.py
@time: 2019/01/08
'''

ruledir="/hwfssz5/ST_BIGDATA/USER/yueyao/snakemake/test/snakemake-example"

include: ruledir+'Rm_rRNA.Filter.snakemake.py'
include: ruledir+'/Filter.circosrna.snakemake_md.py'
