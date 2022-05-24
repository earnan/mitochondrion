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

from sympy import FallingFactorial
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
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='F:\\4228\\pyss\\pyss.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:\\4228\\pyss\\gene.annotation.info', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def name_mapping(s, table=5):  # 名字映射,第5套密码子
    amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                    '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S', 'G']
    amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                    'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Gly']
    #"""cds_dict rrn_dict通用"""
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


def codon_check(codon_str, table=5):  # 第一步,检查cds的起止密码子
    start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
    end_codon_list = ['TAA', 'TAG', 'TA', 'T']
    start_codon = codon_str.split('/')[0]
    end_codon = codon_str.split('/')[-1]
    if start_codon in start_codon_list and end_codon in end_codon_list:
        flag = True
    else:
        flag = False
    return flag


def overlap_check(cds_ovl_dict, gene_list, gene_pos_dict):  # 第二步,检查cds的overlap
    for ovl_cds in cds_ovl_dict.keys():
        ovl_cds_index = gene_list.index(ovl_cds)
        if ovl_cds_index+1 >= len(gene_list):
            next_gene = gene_list[0]
        else:
            next_gene = gene_list[ovl_cds_index+1]
        # print(next_gene)
        if next_gene.startswith('trn'):  # 分四种情况,但只需要考虑cds为正的两种情况
            if gene_pos_dict[ovl_cds].split(':')[-1] == '+':  # 如果cds是正链的话
                start = gene_pos_dict[ovl_cds].split('-')[0]
                new_end = int(gene_pos_dict[next_gene].split('-')[0])-1
                print('{} pos may be {}-{}:+ '.format(ovl_cds, start, new_end))
    return 0


in_path = args.infile
out_path = args.outfile
with open(in_path, 'r') as in_handle, open(out_path, 'w') as out_handle:
    """计数"""
    cds_n = 0
    trn_n = 0
    rrn_n = 0
    dloop_n = 0
    """跳过第一行"""
    line = in_handle.readline()
    tmp_ovl = 0  # 中间变量,用于存储overlap
    cds_ovl_dict = {}  # 用于存储overlap
    gene_list = []  # 存储所有基因
    gene_pos_dict = {}  # 存储所有基因及其位置
    for line in in_handle:
        """分割赋值"""
        line_content = line.split('\t')
        name = line_content[0]
        start = line_content[1]
        end = line_content[2]
        strand = line_content[3]
        lenth = line_content[4]
        ovl = line_content[5]
        codon_str = line_content[6]
        """判断"""
        if line.startswith('rrn'):  # rrna
            rrn_n += 1
            n = 'rRNA'+str(rrn_n)
            name1 = name
            name2 = name_mapping(name)
            gene_list.append(name1)
            gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
            s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                n, start, end, strand, name1, name2)
        elif line.startswith('trn'):  # trna
            trn_n += 1
            n = 'tRNA'+str(trn_n)
            name1 = name.strip(')').split(
                '(')[0]+'-'+name.strip(')').split('(')[1].upper()
            letter = (re.search(r'[A-Z]', name1.split('-')[0])).group(0)
            name2 = 'tRNA-'+name_mapping(letter)
            gene_list.append(name1)
            gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
            s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                n, start, end, strand, name1, name2)
        elif line.startswith('OH'):  # dloop
            dloop_n += 1
            n = 'D-loop'+str(dloop_n)
            gene_list.append(name)
            gene_pos_dict[name] = '{0}-{1}:{2}'.format(start, end, strand)
            s = '{0}\t{1}-{2}:{3}'.format(
                n, start, end, strand)
        else:  # cds
            cds_n += 1
            n = 'CDS'+str(cds_n)
            if int(ovl) < 0:
                cds_ovl_dict[name] = ovl
            name1 = name
            name2 = name_mapping(name.split('-')[0])
            gene_list.append(name1)
            gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
            flag = codon_check(codon_str)
            if flag == True:
                s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                    n, start, end, strand, name1, name2)
            elif flag == False:
                s = '{0}\t{1}-{2}:{3}\t{4}\t{5}\t{6}'.format(
                    n, start, end, strand, name1, name2, flag)
        """写入"""
        print(s)
        out_handle.write(s+'\n')

print(cds_ovl_dict)
print(gene_list)
print(gene_pos_dict)


overlap_check(cds_ovl_dict, gene_list, gene_pos_dict)

"""写入统计文件"""
with open((args.output+os.sep+'log'), 'w') as f_log:
    f_log.write(
        'gene{0}atp6{0}atp8{0}cob{0}cox1{0}cox2{0}cox3{0}nad1{0}nad2{0}nad3{0}nad4{0}nad4L{0}nad5{0}nad6\n'.format('\t'))  # 小写
    f_log.write(
        'gene{0}ATP6{0}ATP8{0}CYTB{0}COX1{0}COX2{0}COX3{0}ND1{0}ND2{0}ND3{0}ND4{0}ND4L{0}ND5{0}ND6\n'.format('\t'))  # 大写
all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                       'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                       'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
"""统计初始化"""
dict_missing_gene = {}  # 每个文件中缺失的基因统计,总 字典
dict_gene_len = {}  # 统计每个基因在不同物种中的长度,取平均
for i in all_gene_list_upper:
    dict_gene_len[i] = []

"""初始化"""
dict_file_cds_count = {}  # 每个文件中cds计数
file_list = os.listdir(args.input)
file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字
if os.path.exists(args.output) == False:
    os.mkdir(args.output)

"""主程序"""
for file in file_list:
    ingbk_path = os.path.join(args.input, file)
    cds_fasta, complete_fasta, count, file_name,  list_gene_name, s, dict_gene_len, seq_id = get_cds(
        ingbk_path, False, dict_gene_len)
    dict_file_cds_count[seq_id] = count  # 每个文件中cds计数
    """写入文件"""
    with open((args.output+os.sep+seq_id+'.fasta'), 'wb') as f_complete, \
            open((args.output+os.sep+file_name.rstrip('.gbk')+'_cds.fasta'), 'wb') as f_cds, \
            open((args.output+os.sep+'log'), 'a+') as f_log:
        f_cds.write(cds_fasta.encode())
        f_complete.write(complete_fasta.encode())

        """以下为统计部分"""
        f_log.write(s+'\n')
        f_log.write('>'+file.rstrip('.gbk')+'\t')
        list_missing_gene = []  # 每个文件中缺失的基因统计,单列表
        for i in range(len(all_gene_list_upper)):
            if all_gene_list_upper[i] in list_gene_name \
                    or all_gene_list_upper[i].lower() in list_gene_name \
                    or all_gene_list_lower[i] in list_gene_name \
                    or all_gene_list_lower[i].upper() in list_gene_name:
                f_log.write(all_gene_list_upper[i]+'\t')
            else:
                f_log.write('NULL'+'\t')
                list_missing_gene.append(all_gene_list_upper[i])  # 缺失的基因
        print(list_missing_gene)
        f_log.write('\n')
        [f_log.write(tmp+'\t') for tmp in list_missing_gene]
        f_log.write('\n')
    dict_missing_gene['>'+seq_id] = list_missing_gene  # 放入字典当中
with open((args.output+os.sep+'log'), 'a+') as f_log:
    f_log.write(str(dict_missing_gene))
