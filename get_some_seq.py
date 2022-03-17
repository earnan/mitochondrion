#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   get_some_seq.py
#         Author:   yujie
#    Description:   get_some_seq.py
#        Version:   1.0
#           Time:   2022/03/17 09:42:02
#  Last Modified:   2022/03/17 09:42:02
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
"""待解决问题
    1.速度较慢.后续查找耗时原因
    2.自定义输出路径
"""
from Bio import SearchIO
import argparse
from Bio import SeqIO
import os
import re
import time
import linecache
import random

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   get_some_seq.py')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='fasta path', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-n', '--number',
                      metavar='[int]', help='抽取比例', type=int, default=3, required=False)

optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

##################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))

#######################################################################################


def read_fasta_to_dic(infasta):
    with open(infasta, 'r') as f:
        seq_id = ''
        seq_index = []
        dict_seq = {}
        dict_len = {}
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip('\n')
                seq_index.append(line.replace(
                    "\n", "").replace(">", ""))
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{2} Item Quantity: {0} {1}'.format(
        len(dict_seq), len(dict_len), os.path.basename(infasta)))
    return dict_seq, dict_len


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note + '\n' + format_seq + "\n"


# 以二进制方式读取文件,统计换行符个数,即行数
count = 0
thefile = open(args.infile, 'rb')
while True:
    buffer = thefile.read(1024 * 8192)
    if not buffer:
        break
    count += buffer.count('\n'.encode())
thefile.close()
print('Item Total {}'.format(count/2))

# 取出1/3的内容写入新文件
dict_all_contents, dict_len = read_fasta_to_dic(args.infile)
num_list = []
dict_content = {}
while True:
    n = 2*random.randint(1, count/2)-1  # >ID
    if n not in num_list:
        num_list.append(n)
        id = linecache.getline(args.infile, n).strip()
        content = linecache.getline(args.infile, n+1).strip()
        dict_content[id] = content
    if len(num_list) > count/2/args.number:
        print('抽取了{}条序列'.format(len(num_list)))
        break
# print(dict_content)
#new = list(set(all_contents)-set(num_content_list))

write_path1 = os.path.join(os.path.dirname(args.infile),
                           'remainder_'+os.path.basename(args.infile))
with open(write_path1, 'w') as f:
    for i in dict_all_contents.keys():
        if i not in dict_content.keys():
            fasta = format_fasta(i, dict_all_contents[i], 70)
            f.write(fasta)
print('剩余序列存放于{}'.format(write_path1))

write_path2 = os.path.join(os.path.dirname(args.infile),
                           'get_'+os.path.basename(args.infile))
with open(write_path2, 'w') as f:
    for i in dict_content.keys():
        fasta = format_fasta(i, dict_content[i], 70)
        f.write(fasta)
print('抽取出序列存放于{}'.format(write_path2))

# print('\n')
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
###########
