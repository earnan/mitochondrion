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
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='F:\\4228\\pyss\\pyss.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:\\4228\\pyss\\gene.annotation.info', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def name_mapping(s, table=5):  # 名字映射,第5套密码子,name2
    amino_acid_1 = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y',
                    '*', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'S', 'G']
    amino_acid_2 = ['Phe', 'Leu', 'Ile', 'Met', 'Val', 'Ser', 'Pro', 'Thr', 'Ala',
                    'Tyr', 'Ter', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Ser', 'Gly']
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
    pw('--------------------------Step 2 Check cds overlap!--------------------------')
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
                pw('{} pos may be {}-{}:+ '.format(ovl_cds, start, new_end))
    return 0


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
    all_trna_list = ['trnK-TTT', 'trnV-TAC', 'trnL1-TAG', 'trnA-TGC', 'trnP-TGG', 'trnD-GTC', 'trnC-GCA', 'trnF-GAA', 'trnY-GTA', 'trnW-TCA',
                     'trnG-TCC', 'trnH-GTG', 'trnQ-TTG', 'trnL2-TAA', 'trnN-GTT', 'trnR-TCG', 'trnE-TTC', 'trnM-CAT', 'trnS2-TGA', 'trnS1-GCT', 'trnT-TGT', 'trnI-GAT']
    all_cds_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                          'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_cds_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                          'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    """初始化"""
    list_missing_cds,  list_missing_trna = [], []  # 缺失
    list_extra_cds, list_extra_trna = [], []  # 多余
    """统计缺失的"""
    for i in range(len(all_cds_list_upper)):
        if not (all_cds_list_upper[i] in cds_list
                or all_cds_list_upper[i].lower() in cds_list
                or all_cds_list_lower[i] in cds_list
                or all_cds_list_lower[i].upper() in cds_list):
            list_missing_cds.append(all_cds_list_upper[i])
    for i in range(len(all_trna_list)):
        if all_trna_list[i] not in trna_list:
            list_missing_trna.append(all_trna_list[i])
    """统计多余的"""
    for i in cds_list:
        if len(i.split('-')) > 1:
            list_extra_cds.append(i)
    for i in trna_list:
        if len(i.split('_')) > 1:
            list_extra_trna.append(i)
    """输出提示"""
    """
    if len(list_missing_cds) > 0:  # 少的cds
        print(list_missing_cds)
    if len(list_missing_trna) > 0:  # 少的trna
        print(list_missing_trna)
    if len(list_extra_cds) > 0:  # 多的cds
        print(list_extra_cds)
    if len(list_extra_trna) > 0:  # 多的trna
        print(list_extra_trna)
    """
    return list_missing_cds, list_missing_trna, list_extra_cds, list_extra_trna


def tbl_format_parse(in_path=args.infile, out_path=args.outfile):
    with open(in_path, 'r') as in_handle, open(out_path, 'w') as out_handle:
        """计数"""
        count, cds_n, trn_n, rrn_n, dloop_n = 0, 0, 0, 0, 0
        """跳过第一行"""
        line = in_handle.readline()
        cds_ovl_dict = {}  # 用于存储有overlap的cds
        gene_list = []  # 存储所有基因
        gene_pos_dict = {}  # 存储所有基因及其位置
        gene_lenth_dict = {}  # 存储所有基因长度的字典
        tmp_line_number_list = []  # 起止密码子有问题的行数
        for line in in_handle:
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
            else:  # cds
                cds_n += 1
                n = 'CDS'+str(cds_n)
                if int(ovl) < 0:
                    cds_ovl_dict[name] = ovl
                name1 = name
                name2 = name_mapping(name.split('-')[0])
                gene_list.append(name1)
                gene_lenth_dict[name1] = lenth
                gene_pos_dict[name1] = '{0}-{1}:{2}'.format(start, end, strand)
                flag = codon_check(codon_str)
                if flag == True:
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}'.format(
                        n, start, end, strand, name1, name2)
                elif flag == False:
                    tmp_line_number_list.append(count)
                    s = '{0}\t{1}-{2}:{3}\t{4}\t{5}\t{6}'.format(
                        n, start, end, strand, name1, name2, flag)
            """写入"""
            print(s)
            out_handle.write(s+'\n')
    return cds_ovl_dict, gene_list, gene_pos_dict, cds_n,    trn_n,    rrn_n,    dloop_n, gene_lenth_dict, count, tmp_line_number_list


def pw(s, out_path=args.outfile):
    print(s)
    with open(out_path, 'a') as out_handle:  # 追加写
        out_handle.write(str(s)+'\n')


print('\n')
cds_ovl_dict, gene_list, gene_pos_dict, cds_n,    trn_n,    rrn_n,    dloop_n, gene_lenth_dict, count, tmp_line_number_list = tbl_format_parse()  # 1
pw('\n')
pw('--------------------------Step 1 Check flag tag!--------------------------')
pw(str(tmp_line_number_list))
overlap_check(cds_ovl_dict, gene_list, gene_pos_dict)  # 2
pw('--------------------------Step 3 Check gene quantity!--------------------------')
list_missing_cds, list_missing_trna, list_extra_cds, list_extra_trna = gene_count_check(
    gene_list)  # 3
pw((cds_n,    trn_n,    rrn_n,    dloop_n))
pw(cds_n+trn_n+rrn_n)
pw(('{}=13+{}-{}'.format(cds_n, len(list_extra_cds), len(list_missing_cds)),
   list_extra_cds, [gene_lenth_dict[i] for i in list_extra_cds], list_missing_cds))
pw(('{}=22+{}-{}'.format(trn_n, len(list_extra_trna),                            len(list_missing_trna)),
   list_extra_trna, [gene_lenth_dict[i] for i in list_extra_trna], list_missing_trna))
pw('--------------------------Step 4 Edit gene name!--------------------------\n--------------------------Step 5 Search rrnl 1&2!--------------------------\n--------------------------Step 6 Search D-loop region!--------------------------')
