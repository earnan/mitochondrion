#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_move_gene_pos_v4.0.py
#         Author:   yujie
#    Description:   mt_move_gene_pos_v4.0.py
#        Version:   1.0
#           Time:   2022/06/07 09:38:21
#  Last Modified:   2022/11/24 17:13:37
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
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
\n\
##########################################################\n\
#\n\
#       Filename:   mt_move_gene_pos.py\n\
#         Author:   yujie\n\
#    Description:   mt_move_gene_pos.py\n\
#        Version:   2.0\n\
#           Time:   2022/06/07 09:38:21\n\
#  Last Modified:   2022/11/24 17:14:32\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   mt_move_gene_pos.py\n\
功能：线粒体基因平移位置/fasta调整起点\n\
1.调整fasta起点,输入 -fa -n2/-s \n\
\tn2:输入正值,将末尾(n2)bp的碱基挪到开头\n\
\ts:输入正值,将开头(s-1)bp的碱基挪到末尾\n\
2.对ann.info2仅排序,输入 -i -o \n\
3.对ann.info2平移+排序\n\
\t不分段操作,输入 -i -o -n2 \n\
\t分段操作,输入 -i -o -ln -n1 -n2 -m\n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\mitochondrion\annotation\mt_move_gene_pos.py\n\
Path: /share/nas1/yuj/script/mitochondrion/annotation/mt_move_gene_pos.py\n\
Version: 2.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--ininfo', metavar='[file]', type=str, help='要修改的注释,使用时需输入', default='', required=False)
optional.add_argument(
    '-fa', '--infasta', metavar='[file]', type=str, help='要修改的fsa,使用时需要输入', default='', required=False)
optional.add_argument(
    '-o', '--outinfo', metavar='[file]', type=str, help='输出文件,使用时需输入', default='', required=False)
optional.add_argument(
    '-ln', '--line_number', metavar='[int]', type=int, help='从第几行开始分段操作,默认1不分段,使用时需输入', default=1,  required=False)
optional.add_argument(
    '-n1', '--number1', metavar='[int]', type=int, help='第一段的平移距离,默认0,使用时需输入', default=0,  required=False)
optional.add_argument(
    '-n2', '--number2', metavar='[int]', type=int, help='第二段的平移距离,默认0,使用时需输入', default=0,  required=False)
optional.add_argument(
    '-s', '--start', metavar='[int]', type=int, help='现有序列的起点,默认0,使用时需输入', default=0,  required=False)
optional.add_argument(
    '-m', '--maxlen', metavar='[len(fasta)]', type=int, help='基因组长度,默认0,使用时需输入', default=0,  required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()

if args.info:
    print(os.path.abspath(__file__))
    print("\n以下是更新日志: ")
    print("\t适用于内含子外显子")
    print("\t20220210\t第三版改为子函数并考虑分段操作")
    print("\t20220607\t增加了排序函数")
    print("\t20220607\t考虑跨首尾基因排序到开头,其自己的计数")
    print("\t20221124\t大幅修改代码运行逻辑,更新使用方法")
    print("\n")
    sys.exit(0)


def multi_replace(origin_string, pos_list, new_str_list):
    '''
    替换多个指定位置上的字符(原字符串,位置列表,新的字符列表),使用时位置列表按自然语言习惯
    把 1-10:- 变为 1-10:minus
    '''
    origin_string_list = []
    for s in origin_string:
        origin_string_list.append(s)
    new_string_list = origin_string_list
    # enumerate(sequence, [start=0]) start定义角标从哪开始,默认为0,同索引计数一致
    for i, pos_info in enumerate(pos_list):
        new_string_list[pos_info-1] = new_str_list[i]
    new_string = ''.join(new_string_list)
    return new_string


def edit_pos(pos_info, n):
    '''
    位置平移
    1-10:+ 修改为 101-110:+
    '''
    if pos_info.endswith('-'):  # 仅针对负链方向,会被修改为minus
        pos_info = multi_replace(pos_info, [0], ['minus'])  # 索引-1对应自然的0
    pos_content = pos_info.split('-')
    new_start = int(pos_content[0])+n  # 例如1
    end_pos_all = pos_content[1]  # 例如10:+
    end_pos_content = end_pos_all.split(':')
    new_end = int(end_pos_content[0])+n  # 例如10
    new_pos_info = (str(new_start)+'-'+str(new_end)+':' +
                    end_pos_content[1]).replace('minus', '-')
    return new_pos_info  # 101-110:+


def get_new_line(line, n, max_len):
    # 20220607  max_len
    line_content = line.split()
    pos_info = line_content[1]  # 例如67574-67687:-;135718-135975:+
    s_content = pos_info.split(';')  # 分号隔开
    if len(s_content) == 1:
        new_pos_info = edit_pos(pos_info, n)
        # 20220610考虑没有内含子的跨首尾的基因
        if max_len != 0 and int(new_pos_info.split('-')[0]) > max_len:
            start = '1'
            end = new_pos_info.split('-')[1]
            if new_pos_info.endswith('-'):  # split会把末尾的 负链也当成分隔符
                new_pos_info = start+'-'+end+'-'
            else:
                new_pos_info = start+'-'+end
    elif len(s_content) == 2:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        # 20220607考虑跨首尾的基因
        if max_len != 0 and int(new_pos_info1.split('-')[0]) > max_len:
            if new_pos_info2.endswith('-'):
                new_pos_info = '1'+'-'+new_pos_info2.split('-')[1]+'-'
            else:
                new_pos_info = '1'+'-'+new_pos_info2.split('-')[1]
        else:
            new_pos_info = new_pos_info1+';'+new_pos_info2
    elif len(s_content) == 3:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        new_pos_info3 = edit_pos(s_content[2], n)
        new_pos_info = new_pos_info1+';'+new_pos_info2+';'+new_pos_info3
    elif len(s_content) > 3:  # 20221124 真菌线粒体基因有很多段
        new_pos_info = ';'.join([edit_pos(s_content[i], n)
                                for i in range(len(s_content))])
    new_line = line.replace(pos_info, new_pos_info)
    print(new_line.strip('\n'))
    return new_line


def pos_sort(input_file, output_file):
    '''
    对注释进行排序并输出
    '''
    with open(input_file, 'r') as fi, open(output_file, 'w') as fo:  # 打开tmp
        dic = {}
        for line in fi:
            id = int(line.split()[1].split('-')[0])  # 位置的第一个起点
            dic[id] = line
        # 20220607考虑跨首尾基因排序到开头,其自己的计数
        count, cds_n, trn_n, rrn_n, dloop_n, ol_n = 0, 0, 0, 0, 0, 0
        """判断"""
        for i in sorted(dic):  # 排序
            line = dic[i].strip('\n')  # 去除了两端换行
            if line.startswith('rRNA'):  # rrna
                rrn_n += 1
                n = 'rRNA'+str(rrn_n)
            elif line.startswith('tRNA'):  # trna
                trn_n += 1
                n = 'tRNA'+str(trn_n)
            elif line.startswith('D-loop'):  # dloop
                dloop_n += 1
                n = 'D-loop'
            elif line.startswith('OL') or line.startswith('rep_origin'):  # ol区  复制起始区域
                ol_n += 1
                n = 'rep_origin'
            elif line.startswith('CDS'):  # cds
                cds_n += 1
                n = 'CDS'+str(cds_n)
            old_str = line.split('\t')[0]
            new_str = n
            print(line.replace(old_str, new_str))
            fo.write(line.replace(old_str, new_str)+'\n')
    return 0


def correcting_fasta(fasta_file, number, start):
    '''
    修改序列起点
    '''
    fasta_path = os.path.abspath(fasta_file)
    indir_path = os.path.dirname(fasta_path)
    file_prefix = os.path.basename(fasta_path).split('.')[0]
    with open(fasta_path, 'r') as fi_handle:
        tmp_dict = {}
        seq_id = fi_handle.readline()
        seq = fi_handle.readline().strip()
        max_len = len(seq)
        tmp_dict[seq_id] = seq
    if number == 0:
        number = max_len-start+1
    """末尾挪到开头"""
    last = seq[-number:]
    # print(last+seq.rstrip(last)) #可能会产生bug
    s = ''
    tmp_list = list(seq)
    for i in range(number):
        tmp_list.pop()  # pop函数默认返回被删除的值  直接用就好
    for i in tmp_list:
        s += i
    new_seq = last+s
    with open(os.path.join(indir_path, file_prefix+'.fsa2'), 'w') as fo_handle:
        fo_handle.write(seq_id)
        fo_handle.write(new_seq+'\n')
    return os.path.join(indir_path, file_prefix+'.fsa2')


if args.infasta != '':  # fasta文件 仅考虑把末尾n2 bp碱基挪到开头
    fsa2_path = correcting_fasta(args.infasta, args.number2, args.start)
    print(fsa2_path)
    print('Done')


if args.ininfo != '' and args.outinfo != '':  # info 排序注释
    if args.number2 > 0 or args.number2 < 0:
        ln = args.line_number
        n1 = args.number1
        n2 = args.number2
        max_len = args.maxlen
        abs_path = os.path.abspath(args.ininfo)
        indir_path = os.path.dirname(abs_path)
        tmp2_outinfo_path = os.path.join(indir_path, 'tmp2')
        tmp_outinfo_path = os.path.join(indir_path, 'tmp')
        with open(args.ininfo, 'r') as fi, open(tmp2_outinfo_path, 'w') as tmp2:
            for i in range(0, ln-1):  # 前(ln-1)行平移n1距离
                line = fi.readline()
                new_line = get_new_line(line, n1, max_len)
                tmp2.write(new_line)
            for line in fi:  # 从第ln行开始平移n2距离
                if line.strip() != '':  # 20220609 考虑输入的注释信息 下面有几行空行
                    new_line = get_new_line(line, n2, max_len)
                    tmp2.write(new_line)
            print('\n')
        pos_sort(tmp2_outinfo_path, tmp_outinfo_path)
    elif args.number2 == 0:  # 直接对已有位置信息的注释进行排序
        #out_dir_path = os.path.dirname(args.outinfo)
        tmp_outinfo_path = args.outinfo+'_tmp'
        pos_sort(args.ininfo, tmp_outinfo_path)
    """
    判断注释文件格式是否有问题 
    判断trna一行是否有问题
    """
    with open(tmp_outinfo_path, 'r') as tmp_handle:
        for line in tmp_handle:
            # 换行符也是一个元素  tRNA1 \t 17-91:- \t tRNA-His \t \n
            # 20220708  重新修改判断,缺少反密码子
            if line.startswith('tRNA') and (not line.split('\t')[2].split('-')[-1].isupper()):
                trn_flag = False
                break
            elif line.startswith('tRNA') and line.split('\t')[2].split('-')[-1].isupper():
                trn_flag = True
                break
    if trn_flag == True:  # 没问题 写进输出文件
        with open(tmp_outinfo_path, 'r') as tmp_handle, open(args.outinfo, 'w') as outinfo_handle:
            content = tmp_handle.read()
            outinfo_handle.write(content)
    elif trn_flag == False:
        print('Please check!')
