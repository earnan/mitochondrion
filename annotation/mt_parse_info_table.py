#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_parse_info_table.py
#         Author:   yujie
#    Description:   mt_parse_info_table.py
#        Version:   1.0
#           Time:   2022/05/23 11:51:11
#  Last Modified:   2022/05/23 11:51:11
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
\npython3   mt_parse_info_table.py\n\
解析网页表格里的注释信息,生成gene.annotation.info\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='F:\\4228\\nd06\\06.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:\\4228\\nd06\\gene.annotation.info', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def name_mapping(s, table=5):  # 第5套密码子
    amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                    '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S', 'G']
    amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                    'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Gly']
    cds_dict = {'nad2': 'NADH dehydrogenase subunit 2',
                'cox1': 'cytochrome c oxidase subunit 1',
                'cox2': 'cytochrome c oxidase subunit 2',
                'atp8': 'ATP synthase F0 subunit 8',
                'atp6': 'ATP synthase F0 subunit 6',
                'cox3': 'cytochrome c oxidase subunit 3',
                'nad3': 'NADH dehydrogenase subunit 3',
                'nad5': 'NADH dehydrogenase subunit 5',
                'nad4': 'NADH dehydrogenase subunit 4',
                'nad4l': 'NADH dehydrogenase subunit 4L',
                'nad6': 'NADH dehydrogenase subunit 6',
                'cob': 'cytochrome b',
                'nad1': 'NADH dehydrogenase subunit 1'}
    rrn_dict = {'rrnL': 'l-rRNA', 'rrnS': 's-rRNA'}
    if len(s) == 1:
        for i in range(len(amino_acid_1)):
            if amino_acid_1[i] == s:
                c = amino_acid_2[i]
    elif s.startswith('rrn'):
        c = rrn_dict[s]
    else:
        c = cds_dict[s]
    return c


in_path = args.infile
out_path = args.outfile
with open(in_path, 'r') as in_handle, open(out_path, 'w') as out_handle:
    cds_n = 0
    trn_n = 0
    rrn_n = 0
    dloop_n = 0
    line = in_handle.readline()
    for line in in_handle:
        line_content = line.split('\t')
        name = line_content[0]
        start = line_content[1]
        end = line_content[2]
        strand = line_content[3]
        lenth = line_content[4]
        nc = line_content[5]
        codons = line_content[6]
        if line.startswith('rrn'):
            rrn_n += 1
            n = 'rRNA'+str(rrn_n)
            name1 = name
            name2 = name_mapping(name)
        elif line.startswith('trn'):
            trn_n += 1
            n = 'tRNA'+str(trn_n)
            name1 = name.strip(')').split(
                '(')[0]+'-'+name.strip(')').split('(')[1].upper()
            letter = (re.search(r'[A-Z]', name1.split('-')[0])).group(0)
            name2 = 'tRNA-'+name_mapping(letter)
        elif line.startswith('OH'):
            dloop_n += 1
            n = 'D-loop'+str(dloop_n)
            name1 = ''
            name2 = ''
        else:
            cds_n += 1
            n = 'CDS'+str(cds_n)
            name1 = name
            name2 = name_mapping(name.split('-')[0])
        s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(n,
                                                start, end, strand, name1, name2)
        print(s)
        out_handle.write(s+'\n')
