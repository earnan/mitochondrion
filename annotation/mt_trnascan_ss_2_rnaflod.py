#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_trnascan_ss_2_rnaflod.py
#         Author:   yujie
#    Description:   mt_trnascan_ss_2_rnaflod.py
#        Version:   1.0
#           Time:   2022/06/06 17:44:54
#  Last Modified:   2022/06/06 17:44:54
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import re
import time
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_trnascan_ss_2_rnaflod.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='F:\\tRNA.ss', required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[outdir]', help='outdir', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

"""
#这段在此程序没用,就是个思路
def read_fasta_to_dic1(infasta):
    with open(infasta, 'r') as f:
        seq_id = ''
        id_index = []
        dict_seq = {}
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip('\n')
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))
                dict_seq[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
    print(len(dict_seq))
    return dict_seq
"""

# ################################################存进一个大字典里
with open(args.infile, 'r') as fi_handle:
    line_n = 0
    multiple = 1
    tmp_dict = {}
    tmp_dict[multiple] = []
    for line in fi_handle:
        line_n += 1
        tmp_dict[multiple].append(line)
        while line_n > 6*multiple:
            multiple += 1
            tmp_dict[multiple] = []

# #########################################存进不同的字典里,需要用到local()函数
