#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_blast_match.py
#         Author:   yujie
#    Description:   mt_blast_match.py
#        Version:   1.0
#           Time:   2022/03/11 09:49:54
#  Last Modified:   2022/03/11 09:49:54
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
from Bio import SeqIO
import os
import re

from sympy import content
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   mt_from_gbk_get_cds.py')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='sample-cds', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-r', '--ref',
                      metavar='[file]', help='ref-cds', type=str, default="F:\\ref_tre\\gene\\feature\\ref.gene.seq", required=False)
optional.add_argument('-o', '--outdir',
                      metavar='[dir]', help='输出的路径', type=str, default="F:\\ref_tre\\gene\\blast", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


cmd = "formatdb -i {} -p F -o F".format(args.infile)
# print(cmd)
# os.system(cmd)

cmd = "blastn -query {0} -db {1} -evalue 1e-5 -outfmt 5 -max_hsps 1 -out {2}/blastn.result.xml -word_size 7".format(
    args.ref, args.infile, args.outdir)
# print(cmd)
# os.system(cmd)

cmd = "perl /share/nas6/xul/program/mt2/phytree/gene_tree/src/blast_parser.pl -tophit 1 -topmatch 1 -m 7 {0}/blastn.result.xml > {0}/blastn.tophit.result.xls".format(
    args.outdir)
# print(cmd)
# os.system(cmd)


"""
根据blast 结果
输出infile每一行信息
输出ref对应长度信息
把提及到的序列放一起
"""
blastn_tophit_result_path = os.path.join(
    args.outdir, 'blastn.tophit.result.xls')
#fi = open(blastn_tophit_result_path, 'r')
with open(blastn_tophit_result_path, 'r') as f:
    for i in range(0, 259):
        f.readline()
    for line in f:
        print(line.strip())
        content = line.split('\t')

        print('>{0} len: {1}'.format(content[-1].strip(), content[5].strip()))
        #subject_start_end_list = re.findall(r'\d+', content[15].split()[1])

        print('{0} {1}'.format(content[0].split()[0], content[1].strip()))
        #query_start_end_list = re.findall(r'\d+', content[0].split()[1])
