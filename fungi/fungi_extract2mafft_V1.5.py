#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   fungi_extract2mafft_V1.5.py
#         Author:   yujie
#    Description:   fungi_extract2mafft_V1.5.py
#        Version:   2.0
#           Time:   2023/02/23 09:23:02
#  Last Modified:   2023/02/23 09:23:02
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
import copy  # 深度拷贝
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
Fungi:\n\
    step1:python3 fungi_extract2mafft_V1.5.py -i cds_all/ -o1 gene/extract\n\
    step2:python3 fungi_extract2mafft_V1.5.py -i cds_all/ -o1 gene/extract -o2 gene/mafft -c\n\
Usually:\n\
    python3 fungi_extract2mafft_V1.5.py -i cds_all/ -o1 gene/extract -o2 gene/mafft\n\
    提取同名基因序列\n\
    mafft比对\n\
    v2.0并行运行')
optional = parser.add_argument_group('optional')
required = parser.add_argument_group('required')
required.add_argument('-i', '--input',
                      metavar='[dir]', help='cds file dir path', type=str)
optional.add_argument('-o1', '--outdir1',
                      metavar='[dir]', help='gene/extract dir path', type=str)
optional.add_argument('-o2', '--outdir2',
                      metavar='[dir]', help='gene/mafft dir path', type=str)
optional.add_argument('-c', '--check', help='check target gene,default False',
                      action='store_true', required=False)
optional.add_argument('-ignore', '--ignore', help='cancel ignore warning message,default true',
                      action='store_false', required=False)
optional.add_argument('-info', '--info', help='show update log and exit',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20220906 修改基因组的类型判断,若不符合则退出')
    print('\t20221012 修改基因名字映射函数,使其能够识别CO1样式')
    print('\t20221111 增加gbk里基因名是全称(如ATPase)时的处理')
    print('\t20221217 feat: ✨ 对gbk文件进行去重')
    print('\t20221219 🐞fix(get_gene_note): 修改没有gene标签的cds的默认ID')
    print('\t20221220 ✨feat(main): 写入去重后的登录号')
    print('\n')
    sys.exit(0)


#################################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
print(f'Start Time : {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}\n')
#################################################################


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        # format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


def read_fasta_to_dic1(infasta):  # 最简单,针对普通fasta文件 >物种名
    with open(infasta, 'r') as f:
        seq_id = ''
        id_index = []
        dict_seq = {}
        for line in f:
            # 如果是">"开头的，就创建一个key键
            if line.startswith('>'):
                seq_id = line.strip('\n')  # ID为键
                id_index.append(line.replace(
                    "\n", "").replace(">", ""))  # 顺便创建索引的列表
                dict_seq[seq_id] = ''  # 有key无value
            # 如果不是">"开头的，在这个key键下添加value值
            else:
                dict_seq[seq_id] += line.strip('\n')  # 包含说明行序列行的字典结构
    length_dict_seq = len(dict_seq)
    # print(length_dict_seq)
    return dict_seq, length_dict_seq


# ##################################################格式化基因名字,可重复使用,首次出现于mt_from_gbk_get_cds.py
def gene_name_standardization(gene_name):
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'RPS3', 'ATP9']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6', 'rps3', 'atp9']
    gene_name = gene_name.replace('III', '3').replace(
        'II', '2').replace('I', '1')  # 20220825 gbk文件里 基因名是罗马数字，因此需要先处理 #20220831 copy过来
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 15:
            if all_gene_list_lower[i] == gene_name.lower():
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 15:
            if args.ignore == False:
                print(gene_name, 'WARNING!Please check!')
            gene_name = False
    return gene_name


# ###################################################主函数
if __name__ == '__main__':

    if (args.outdir1) and (args.check == False):
        if not os.path.exists(args.outdir1):
            os.makedirs(args.outdir1)
        all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                               'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'RPS3', 'ATP9']
        all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                               'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6', 'rps3', 'atp9']
        dict_gene_id_seq = {}
        for i in all_gene_list_upper:
            dict_gene_id_seq[i] = ''

        file_list = os.listdir(args.input)
        file_list.sort()
        total_file_number = 0
        total_cds_number = 0
        for file in file_list:
            total_file_number += 1
            species_id = file.lstrip('cds_').replace(
                '.fasta', '')  # 20220627 随上一个脚本更改
            dict_seq, length_dict_seq = read_fasta_to_dic1(
                os.path.join(args.input, file))
            s = f'--------------------------------------------------\
--------------------------------------------------\n{total_file_number}:{species_id}:{length_dict_seq}\n'
            total_cds_number += length_dict_seq
            gene_name_dict = {}
            for j in dict_seq.keys():
                gene = j.split()[-1].split('=')[-1].rstrip(']')
                gene = gene_name_standardization(gene)
                if gene != False:
                    if gene not in gene_name_dict.keys():
                        gene_name_dict[gene] = [len(dict_seq[j])]
                    else:
                        gene_name_dict[gene].append(len(dict_seq[j]))
                    dict_gene_id_seq[gene] += format_fasta(
                        j, dict_seq[j], 70)
            display_flag = False
            for k, v in gene_name_dict.items():
                if len(v) > 1:
                    display_flag = True
                    s += f'{str(k).ljust(20)} {str(v)}\n'
            if display_flag == True:
                print(s.strip())
        average_number = total_cds_number/total_file_number

        if average_number > 13:
            all_gene_list_upper2 = all_gene_list_upper
        else:
            all_gene_list_upper2 = all_gene_list_upper[:13]
        n = 0
        for i in all_gene_list_upper2:
            n += 1
            filename = 'gene{0}.{1}.fasta'.format(n, i)
            with open(os.path.join(args.outdir1,  filename), 'wb') as f:
                f.write(dict_gene_id_seq[i].encode())

    if args.outdir2:
        if not os.path.exists(args.outdir2):
            os.makedirs(args.outdir2)
        if os.path.exists('mafft.sh'):
            os.system('rm mafft.sh')
        if os.path.exists('fasta2line.sh'):
            os.system('rm fasta2line.sh')
        #########################################
        file_list1 = os.listdir(args.outdir1)
        file_list1.sort()
        for file1 in file_list1:
            infasta1 = os.path.join(args.outdir1, file1)
            os.system(
                f'echo "mafft --auto --quiet {infasta1} > {args.outdir2}/{file1}.aln"  >> mafft.sh')
    #########################################
        for file1 in file_list1:
            file2 = f'{file1}.aln'
            inaln2 = os.path.join(args.outdir2, file2)
            cmd = f'echo "perl /share/nas6/xul/program/mt2/phytree/gene_tree/src/fasta2line.pl -i {inaln2} \
-o {args.outdir2}/{file2.replace(".aln", "")}" >> fasta2line.sh'
            os.system(cmd)
    #########################################
        print(
            f'Mafft Start : {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}')
        os.system('sh mafft.sh')
        print(
            f'Convert Start : {time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())}')
        os.system('sh fasta2line.sh')
        os.system(f'rm {args.outdir2}/*.aln')

###############################################################
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('--------------------------------------------------\
--------------------------------------------------\n\nEnd Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################
