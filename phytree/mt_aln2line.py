#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_aln2line.py
#         Author:   yujie
#    Description:   mt_aln2line.py
#        Version:   1.0
#           Time:   2022/06/20 17:15:42
#  Last Modified:   2022/06/20 17:15:42
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
\npython3   mt_aln2line.py多行变单行\n\
E:\OneDrive\jshy信息部\Script\mitochondrion\phytree\mt_aln2line.py\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default="F:\\3816\\allpopulation.fa.aln", required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:\\3816\\allpopulation.fasta', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


with open(args.infile, 'r') as fi_handle, open(args.outfile, 'wb') as fo_handle:
    fi_dict = {}
    for line in fi_handle:
        if line.startswith('>'):
            seq_id = line.strip()
            fi_dict[seq_id] = ''
        else:
            fi_dict[seq_id] += line.strip()

    for i in sorted(fi_dict):  # 排序,默认按键排序
        fo_handle.write((i+'\n').encode())
        fo_handle.write((fi_dict[i]+'\n').encode())
