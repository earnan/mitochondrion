#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_parse_info_table.py
#         Author:   yujie
#    Description:   mt_parse_info_table.py
#        Version:   1.0
#           Time:   2022/05/23 11:51:11
#  Last Modified:   2022/11/08 11:12:42
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
# import linecache  # 大文件行读取
# import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
#import time
# import copy  # 深度拷贝
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   mt_parse_info_table.py\n\
#         Author:   yujie\n\
#    Description:   mt_parse_info_table.py\n\
#        Version:   1.0\n\
#           Time:   2022/05/23 11:51:11\n\
#  Last Modified:   2022/11/08 11:14:35\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   mt_parse_info_table.py\n\
解析网页表格里的注释信息,生成gene.annotation.info2\n\
Function:\n\
1.常规使用\n\
1.1 -i [gene.info] -o [gene.annotation.info2] -n [2/5]\n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\mitochondrion\annotation\mt_parse_info_table.py\n\
Path: /share/nas1/yuj/script/mitochondrion/annotation/mt_parse_info_table.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', type=str, default='F:\\4313\\nano\\nano.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', type=str, default='F:\\4313\\nano\\gene.annotation.info', required=False)
optional.add_argument(
    '-n', '--tablenumber', metavar='[codon table default2]', type=int, default=2, required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20221108 vim生成的gene.info文件最后一行为空，要跳出')
    print('\n')
    sys.exit(0)


# ##################################################################################名字映射


def name_mapping(s, table=5):  # 名字映射,第5套密码子,name2  其实 2  5 一样的
    if table == 5:
        amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                        '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S', 'G']
        amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                        'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Gly']
    elif table == 2:
        amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                        '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S', '*', 'G']
        amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                        'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Ter', 'Gly']

    # """cds_dict rrn_dict通用"""
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
        c = cds_dict[s.split('_')[0]]  # 最开始 只考虑 atp-a 形式  后面考虑 nad5_0 形式
    return c


def trna_mapping(name1):  # trna 映射  临时加的子函数   以后有空再写全
    # 不如用 re.findall(r'[A-Z]', name1.split('-')[0])[0]
    letter = (re.search(r'[A-Z]', name1.split('-')[0])).group(0)
    name2 = 'tRNA-'+name_mapping(letter)
    return name2

# ########################################################################################################################step 1


def codon_check(codon_str, table):  # 第一步,检查cds的起止密码子
    # 20220525   考虑	ATG/T(AA)形式
    # 20220601   考虑    2 5起止密码子 不同
    if table == 5:
        start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 5
    elif table == 2:
        start_codon_list = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'AGA',
                          'AGG', 'TA', 'T', 'AG']  # 2,转录时要加A
    start_codon = codon_str.split('/')[0]
    end_codon = codon_str.split('/')[-1]

    if start_codon in start_codon_list:
        if end_codon in end_codon_list:
            flag = 0  # ok
        elif end_codon not in end_codon_list:
            if end_codon.split('(')[0] in end_codon_list:
                flag = 0
            elif end_codon.split('(')[0] not in end_codon_list:
                flag = 1  # 仅末尾x
    elif start_codon not in start_codon_list:
        if end_codon in end_codon_list:
            flag = 2  # 仅开头x
        elif end_codon not in end_codon_list:
            if end_codon.split('(')[0] in end_codon_list:
                flag = 2
            elif end_codon.split('(')[0] not in end_codon_list:
                flag = 3  # 都错

                # if start_codon in start_codon_list and end_codon in end_codon_list:
                #flag = True
                # else:
                # if end_codon.split('(')[0] in end_codon_list:
                #flag = True
                # else:
                #flag = False
    return flag

# ##################################################################################################step 2


def overlap_check(cds_ovl_dict, trn_ovl_dict, gene_list, gene_pos_dict):  # 第二步,检查cds的overlap

    for ovl_cds in cds_ovl_dict.keys():  # cds_ovl_dict
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
                pw('{} pos may be {}-{}:+ '.format(ovl_cds, start, new_end))

    for ovl_trn in trn_ovl_dict.keys():  # 20220727 trn_ovl_dict
        ovl_trn_index = gene_list.index(ovl_trn)
        if ovl_trn_index+1 >= len(gene_list):
            next_gene = gene_list[0]  # 下一个基因是开头
        else:
            next_gene = gene_list[ovl_trn_index+1]
        # 分四种情况,但只需要考虑trna为负的两种情况
        if not next_gene.startswith('OL') and not next_gene.startswith('OH') and not next_gene.startswith('trn') and not next_gene.startswith('rrn'):
            if gene_pos_dict[ovl_trn].split(':')[-1] == '-':  # 如果trna是负链的话
                new_start = int(gene_pos_dict[ovl_trn].split(
                    '-')[-2].split(':')[0])+1
                end = gene_pos_dict[next_gene].split('-')[-2]
                pw('{} pos may be {}-{}- '.format(next_gene, new_start, end))
    return 0

# #########################################################################################step 3


def gene_count_check(gene_list):  # 第三步,检查缺失和多余的基因
    """分成4个列表"""
    cds_list, rrna_list, trna_list, dloop_list = [], [], [], []  # name1 形式 dloop OH_形式
    for i in gene_list:
        if i.startswith('rrn'):
            rrna_list.append(i)
        elif i.startswith('trn'):
            trna_list.append(i)
        elif i.startswith('OH'):
            dloop_list.append(i)
        else:
            cds_list.append(i)
    """赋值4个列表"""
    all_rrna_list = ['rrnL', 'rrnS']
    all_trna_list = ['trnK-TTT', 'trnV-TAC', 'trnL1-TAG', 'trnA-TGC', 'trnP-TGG', 'trnD-GTC', 'trnC-GCA', 'trnF-GAA', 'trnY-GTA', 'trnW-TCA',
                     'trnG-TCC', 'trnH-GTG', 'trnQ-TTG', 'trnL2-TAA', 'trnN-GTT', 'trnR-TCG', 'trnE-TTC', 'trnM-CAT', 'trnS2-TGA', 'trnS1-GCT', 'trnT-TGT', 'trnI-GAT']
    all_cds_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                          'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_cds_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                          'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    """初始化"""
    list_missing_cds,  list_missing_trna, list_missing_rrna = [], [], []  # 缺失
    list_extra_cds, list_extra_trna, list_extra_rrna = [], [], []  # 多余
    """统计缺失的"""
    for i in range(len(all_cds_list_upper)):  # cds
        if not (all_cds_list_upper[i] in cds_list
                or all_cds_list_upper[i].lower() in cds_list
                or all_cds_list_lower[i] in cds_list
                or all_cds_list_lower[i].upper() in cds_list):
            list_missing_cds.append(all_cds_list_upper[i])
    for i in range(len(all_trna_list)):  # trna
        if all_trna_list[i] not in trna_list:
            list_missing_trna.append(all_trna_list[i])
    for i in range(len(all_rrna_list)):  # rrna
        if all_rrna_list[i] not in rrna_list:
            list_missing_rrna.append(all_rrna_list[i])
    """统计多余的"""
    for i in cds_list:
        if len(i.split('-')) > 1:  # 仅考虑atp-a形式
            list_extra_cds.append(i)
        elif len(i.split('_')) > 1:  # 0602 考虑nad5_0形式
            list_extra_cds.append(i)

    for i in trna_list:
        if len(i.split('_')) > 1:
            list_extra_trna.append(i)
    return list_missing_cds, list_missing_trna, list_missing_rrna, list_extra_cds, list_extra_trna

# ########################################################################################################################################


def tbl_format_parse(in_path=args.infile, out_path=args.outfile, table=args.tablenumber):
    with open(in_path, 'r') as in_handle, open(out_path, 'w') as out_handle:
        """初始化"""
        cds_ovl_dict = {}  # 用于存储有overlap的cds
        trn_ovl_dict = {}  # 20220727 用于存储有overlap的trna
        gene_list = []  # 存储所有基因
        gene_pos_dict = {}  # 存储所有基因及其位置
        gene_lenth_dict = {}  # 存储所有基因长度的字典
        tmp_line_number_list = []  # 起止密码子有问题的行数
        """计数"""
        count, cds_n, trn_n, rrn_n, dloop_n, ol_n = 0, 0, 0, 0, 0, 0
        """跳过第一行 需要提示"""
        line = in_handle.readline()
        if not line.startswith('Name'):
            print(
                'Wrong format! Please add head!\nName	Start	Stop	Strand	Length	ovl/nc	Codons	Infos')
            # 终止程序
            sys.exit(0)  # 0：正常退出 1：异常退出
        for line in in_handle:
            if line.strip() == '':  # 20221108 vim生成的gene.info文件最后一行为空，要跳出
                break
            count += 1
            """分割赋值"""
            line_content = line.split('\t')
            name = line_content[0]
            start = line_content[1]
            end = line_content[2]
            strand = line_content[3]
            lenth = line_content[4]  # 基因长度
            ovl = line_content[5]
            codon_str = line_content[6]
            """判断"""
            if line.startswith('rrn'):  # rrna
                rrn_n += 1
                n = 'rRNA'+str(rrn_n)
                name1 = name
                name2 = name_mapping(name)
                gene_list.append(name1)
                gene_lenth_dict[name1] = lenth
                gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
                s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                    n, start, end, strand, name1, name2)
            elif line.startswith('trn'):  # trna
                trn_n += 1
                n = 'tRNA'+str(trn_n)
                name1 = name.strip(')').split(
                    '(')[0]+'-'+name.strip(')').split('(')[1].upper()
                if int(ovl) < 0:
                    trn_ovl_dict[name1] = ovl
                letter = (re.search(r'[A-Z]', name1.split('-')[0])).group(0)
                name2 = 'tRNA-'+name_mapping(letter)
                gene_list.append(name1)
                gene_lenth_dict[name1] = lenth
                gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
                s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                    n, start, end, strand, name1, name2)
            elif line.startswith('OH'):  # dloop
                dloop_n += 1
                n = 'D-loop'+str(dloop_n)
                gene_list.append(name)  # name=OH_1
                gene_lenth_dict[name] = lenth
                gene_pos_dict[name] = '{0}-{1}:{2}'.format(start, end, strand)
                s = '{0}\t{1}-{2}:{3}'.format(
                    n, start, end, strand)
            elif line.startswith('OL') or line.startswith('rep_origin'):  # ol区  复制起始区域
                ol_n += 1
                n = 'OL'+str(ol_n)
                gene_list.append(name)
                gene_lenth_dict[name] = lenth
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
                gene_lenth_dict[name1] = lenth
                gene_pos_dict[name1] = '{0}-{1}:{2}'.format(
                    start, end, strand)
                flag = codon_check(codon_str, table)

                if flag == 0:
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                        n, start, end, strand, name1, name2)
                elif flag == 1:
                    tmp_line_number_list.append(count)
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}\t{6}'.format(
                        n, start, end, strand, name1, name2, 'end wrong')
                elif flag == 2:
                    tmp_line_number_list.append(count)
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}\t{6}'.format(
                        n, start, end, strand, name1, name2, 'start wrong')
                elif flag == 3:
                    tmp_line_number_list.append(count)
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}\t{6}'.format(
                        n, start, end, strand, name1, name2, '!!!wrong!!!')
            """写入"""
            print(s)
            out_handle.write(s+'\n')
    return cds_ovl_dict, trn_ovl_dict, gene_list, gene_pos_dict, cds_n,    trn_n,    rrn_n,    dloop_n, ol_n, gene_lenth_dict, count, tmp_line_number_list


def pw(s, out_path=args.outfile):
    print(s)
    with open(out_path, 'a') as out_handle:  # 追加写
        out_handle.write(str(s)+'\n')


# ###########################################################################################################################################主函数
print('\n')
cds_ovl_dict, trn_ovl_dict, gene_list, gene_pos_dict, cds_n,    trn_n,    rrn_n,    dloop_n, ol_n, gene_lenth_dict, count, tmp_line_number_list = tbl_format_parse()  # 1
pw('\n')
pw('--------------------------Step 1 Check flag tag!--------------------------')
pw(str(tmp_line_number_list))
pw('--------------------------Step 2 Check cds overlap!--------------------------')
overlap_check(cds_ovl_dict, trn_ovl_dict, gene_list, gene_pos_dict)  # 2
pw('--------------------------Step 3 Check gene quantity!--------------------------')
list_missing_cds, list_missing_trna, list_missing_rrna, list_extra_cds, list_extra_trna = gene_count_check(
    gene_list)  # 3

pw('Total:{} CDS:{} tRNA:{} rRNA:{} D-loop:{} OL:{}'.format(cds_n +
   trn_n+rrn_n, cds_n, trn_n,    rrn_n,    dloop_n, ol_n))
pw('CDS:{}=13+{}-{} | Extra:{} Len:{} Missing:{}'.format(cds_n, len(list_extra_cds), len(list_missing_cds),
   list_extra_cds, [gene_lenth_dict[i] for i in list_extra_cds], list_missing_cds))
pw('tRNA:{}=22+{}-{} | Extra:{} Len:{} Missing:{} Anticodon:{}'.format(trn_n, len(list_extra_trna), len(list_missing_trna), list_extra_trna,
   [gene_lenth_dict[i] for i in list_extra_trna], list_missing_trna, [trna_mapping(i) for i in list_missing_trna]))
pw('rRNA:{}=2-{} | Missing:{}'.format(rrn_n,
   len(list_missing_rrna), list_missing_rrna))
pw('--------------------------Step 4 Edit gene name!--------------------------\n--------------------------Step 5 Search rrnl 1&2!--------------------------\n--------------------------Step 6 Search D-loop region!--------------------------')
