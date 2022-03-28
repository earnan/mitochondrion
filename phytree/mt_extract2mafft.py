#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_extract2mafft.py
#         Author:   yujie
#    Description:   mt_extract2mafft.py
#        Version:   1.0
#           Time:   2022/03/28 17:00:58
#  Last Modified:   2022/03/28 17:00:58
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SearchIO
import argparse
from Bio import SeqIO
import os
import re
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_extract2mafft.py\n\
提取同名基因序列\n\
mafft比对')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入cds所在目录', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\out\\cds', required=False)
optional.add_argument('-o1', '--output1',
                      metavar='[dir]', help='序列提取后存放位置', type=str, default='F:\\ref_tre\\gene\\blast\\fasta', required=False)
optional.add_argument('-o2', '--output2',
                      metavar='[dir]', help='比对好的序列', type=str, default='F:\\ref_tre\\gene\\mafft', required=False)
optional.add_argument('-c', '--check',
                      metavar='[bool]', help='是否用GAP构造序列,默认否,使用时-c 1', type=bool, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

#################################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))
#################################################################

"""模板函数"""


def read_fasta_to_dic3(infasta):  # 适用于带详细位置cds的fa文件
    with open(infasta, 'r') as f:
        seq_id = ''  # 基因名
        dict_seq = {}  # 基因名-序列
        dict_len = {}  # 基因名-长度
        dict_pos = {}  # 基因名-位置,sort()
        d_pos = {}  # 基因名-位置,未排序
        for line in f:
            list_gene_pos = []  # 某基因对应的位置组成的列表 sort
            l_gene_pos = []
            list_n = [0]  # 计算同名基因是第几个
            if line.startswith('>'):
                seq_id = line.strip('\n').split()[2].split('=')[
                    1].strip(']')  # 基因名
                if seq_id in dict_seq.keys():
                    list_n.append(0)
                    seq_id = seq_id+'-'+str(len(list_n))  # 基因名+1 ycf1-2形式
                list_tmp = re.findall(
                    r'\d+', line.strip('\n').split()[1].lstrip(
                        '[').rstrip(']'))  # 位置打散成一个个起点或终点
                [l_gene_pos.append(int(i))
                 for i in list_tmp]  # 转换成数字,放进l_gene_pos,未排序
                d_pos[seq_id] = l_gene_pos

                l_gene_pos.sort()
                [list_gene_pos.append(i) for i in l_gene_pos]
                dict_pos[seq_id] = list_gene_pos
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{0} Item Total: {1} {2} {3}'.format(os.path.basename(infasta),
                                               len(dict_seq), len(dict_len), len(dict_pos)))
    return dict_seq, dict_len, dict_pos, d_pos


###############################################################
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################
