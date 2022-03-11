#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_generate_cmd.py
#         Author:   yujie
#    Description:   mt_generate_cmd.py
#        Version:   1.0
#           Time:   2022/03/10 13:49:04
#  Last Modified:   2022/03/10 13:49:04
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
from Bio import SeqIO
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   快速生成查找单基因命令')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[file]', help='输入log所在路径', type=str, default="F:\Hibiscus_sabdariffa\out\log", required=False)
optional.add_argument('-o', '--output',
                      metavar='[dir]', help='输出的路径', type=str, default="F:/Hibiscus_sabdariffa/out", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


file_input = open(args.input, 'r')
file_output = open(os.path.join(args.output, 'find_cds_cmd'), 'w')
for line in file_input:
    if line.startswith('{'):
        tmp_dict = eval(line)
        print(tmp_dict)


for i in tmp_dict.keys():
    print(len(tmp_dict[i]))
    for j in tmp_dict[i]:
        print(j)
        file_output.write(
            'cp_annotation_one_gene_by_ref_gbk2.pl -i1 fasta/{0}.fasta -i2 gbk/NC_014809.1.gbk -g {1}'.format(i.lstrip('>'), j))
        file_output.write('\n')

file_input.close()
file_output.close()
