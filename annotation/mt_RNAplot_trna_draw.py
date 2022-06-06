#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_RNAplot_trna_draw.py
#         Author:   yujie
#    Description:   mt_RNAplot_trna_draw.py
#        Version:   1.0
#           Time:   2022/06/06 17:20:08
#  Last Modified:   2022/06/06 17:20:08
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

from Bio import SeqIO
from Bio.Seq import Seq
#from icecream import ic
import argparse
import linecache
import os
import re
import time
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_RNAplot_trna_draw.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i1', '--infile1', metavar='[file]', help='info文件', type=str, default='F:\\4228\\hpl01\\gene.annotation.info', required=False)
optional.add_argument(
    '-i2', '--infile2', metavar='[file]', help='fasta文件', type=str, default='F:\\4228\\hpl01\\HPL_01_FULLMT.fsa', required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[dir]', help='out目录', type=str, default='F:\\4228\\hpl01\\trna.structure', required=False)

optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


# ########################################################################################读取gene.info文件
def tbl_format_parse2(in_path=args.infile1):
    # , open(out_path, 'w') as out_handle:
    with open(in_path, 'r') as in_handle:
        """计数"""
        count, cds_n, trn_n, rrn_n, dloop_n = 0, 0, 0, 0, 0
        cds_ovl_dict = {}  # 用于存储有overlap的cds
        gene_list = []  # 存储所有基因
        gene_pos_dict = {}  # 存储所有基因及其位置
        gene_lenth_dict = {}  # 存储所有基因长度的字典
        tmp_line_number_list = []  # 起止密码子有问题的行数

        for line in in_handle:
            count += 1
            """判断"""
            if line.startswith('tRNA'):  # trna
                trn_n += 1
                """分割赋值"""
                line_content = line.split('\t')
                gene_type = line_content[0]
                gene_pos = line_content[1]
                gene_name = line_content[2]
                gene_name1 = line_content[3]

                gene_list.append(gene_name)
                gene_pos_dict[gene_name] = gene_pos
    print('--------------------------------------------------info done!----------------------------------------------------------------------')
    return gene_list, gene_pos_dict

# #######################################################################################读取fasta进行下一步


def read_file(infasta):  # 读取文件
    with open(infasta, 'r') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.strip('\n')
    return seq


def format_pos(pos_str):  # 读取输入的位置为位置列表
    pos_list = []
    content = pos_str.split(';')
    for ele in content:
        if ele.split(':')[-1] == '-':
            tmp = ele.split(':')[0]+':'+'-1'
            pos_list.append(tmp)
        elif ele.split(':')[-1] == '+':
            tmp = ele.split(':')[0]+':'+'1'
            pos_list.append(tmp)
    # print(pos_list)
    return pos_list


def ir(s):  # 反向互补
    re = s[::-1]  # 字符串反向
    c = ""  # 定义字符串c接收互补序列
    for i in re:
        if i == 'A':
            c = c + 'T'
        elif i == 'G':
            c = c + 'C'
        elif i == 'T':
            c = c + 'A'
        elif i == 'C':
            c = c + 'G'
    return c


def merge_sequence(pos_list, seq):  # 合并获取到的序列,顺便排一下位置顺序
    """pos_list 某基因格式化的位置"""
    """seq 全长序列"""
    # -----------------20220523 解决跨首尾基因
    seq_len = len(seq)
    if int(pos_list[0].split(':')[-1]) == 1 and int(pos_list[0].split(':')[0].split('-')[0]) > int(pos_list[0].split(':')[0].split('-')[-1]):  # 14323-1527:1
        pos1 = '{0}-{1}:1'.format(pos_list[0].split(':')
                                  [0].split('-')[0], seq_len)
        pos2 = '1-{}:1'.format(pos_list[0].split(':')[0].split('-')[-1])
        pos_list = [pos1, pos2]
    # ------------------
    cds_seq = ""
    if int(pos_list[0].split(':')[-1]) == -1:  # 一般来说,有内含子的基因,几段方向相同,因此只判断第一段,完了重新排序即可
        pos_list = pos_list[::-1]

    for ele in pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        start_index = start-1
        end_index = end
        if strand == (-1):
            # ic('minus')
            # seq[start_index:end_index] 角标从start_index到end_index    取的是索引start-1一直到end  取的是start一直到end的碱基
            cds_seq += ir(seq[start_index:end_index])
            # ic(cds_seq)
        elif strand == (1):
            # ic('plus')
            cds_seq += seq[start_index:end_index]
            # ic(cds_seq)
    # print(pos_list)
    return cds_seq, pos_list


def get_trna_seq(gene_pos_dict, in_path=args.infile2, outdir_path=args.outdir):
    print('----------------------------------------------------------get seq-----------------------------------------------------------')
    seq = read_file(in_path)
    for i in gene_pos_dict.keys():
        pos_list = format_pos(gene_pos_dict[i])
        cds_seq, pos_list = merge_sequence(pos_list, seq)

        tmp_list = pos_list[0].split(':')[0].split('-')
        prefix = 'ss-{}-{}-{}'.format(i.split('-')
                                      [0], tmp_list[0], tmp_list[1])

        if not os.path.exists(os.path.join(outdir_path, 'trn')):
            os.makedirs(os.path.join(outdir_path, 'trn'))  # 多层创建目录
        outfa_path = os.path.join(outdir_path, 'trn', prefix+'.fa')
        outfold_path = os.path.join(outdir_path, 'trn', prefix+'.fold')
        with open(outfa_path, 'wb') as fo:  # 二进制方式
            fo.write('>{}\n{}\n'.format(prefix, cds_seq).encode())
        #print('----------------------------------------------run cmd------------------------------------------------------------------------')
        cmd1 = 'RNAfold  < {}  > {}'.format(outfa_path, outfold_path)
        os.system(cmd1)

        """
        图片改名
        添加trna-反密码子(对该氨基酸所对应密码子的反向互补
        /share/nas1/yuj/project/GP-20211206-3816/archive/analysis/annotation/Mm_G1/trna.structure/final_tRNA
        cd trna.structure/ && nohup perl /share/nas6/xul/program/mt/tRNA/draw_tRNA.pl -i trn*/*.svg
        
        svg2xxx -t pdf trnH.svg
        """
    print('-----------------------------------------------------cmd done---------------------------------------------------------------------------')
    return 0


# #######################################################################################################主函数
gene_list, gene_pos_dict = tbl_format_parse2()
get_trna_seq(gene_pos_dict)

s = "cd trna.structure/trn && for i in *.fold;do echo $i;RNAplot -o svg < $i;done && rename _ss.svg .svg *.svg && rm *.fa *.fold && cd ../../ && rm *.ps"
print(s)
