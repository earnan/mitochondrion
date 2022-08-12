#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_search_dloop_pos.py
#         Author:   yujie
#    Description:   mt_search_dloop_pos.py
#        Version:   1.0
#           Time:   2022/05/07 15:30:37
#  Last Modified:   2022/05/07 15:30:37
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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_search_dloop_pos.py\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='K125.coverage.txt', type=str, default='F:/4228/nd/K125.coverage.txt', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_file(infile):
    dict_high = {}
    dict_low = {}
    dict_all = {}
    test_content = pd.read_table(infile, header=None)  # 没有标题
    list_depth = test_content[3]
    seq_len = len(list_depth)
    sum_n = sum(list_depth)
    ave_n = int(sum_n/seq_len)
    for i, ele in enumerate(list_depth):  # i为索引,coverag文件的第3列也可看成索引,则实际上对应第i+1个位置
        dict_all[i] = ele
        if ele > 2*ave_n:  # 阈值设置成2
            dict_high[i] = ele
        else:
            dict_low[i] = ele
    return seq_len, sum_n, ave_n, dict_high, dict_low, dict_all


if __name__ == '__main__':
    seq_len, sum_n, ave_n, dict_high, dict_low, dict_all = read_file(
        args.infile)
    # ic(dict_high)
    print('seq_len:{}'.format(seq_len))
    print('sum_n:{}'.format(sum_n))
    print('ave_n:{}'.format(ave_n))
    print('--------------------------------------------------')
    # 把所有位置存为列表
    list_all_keys = []
    list_all_values = []
    for k, v in dict_all.items():
        list_all_keys.append(k)
        list_all_values.append(v)
    # 把覆盖深度高的位置存为列表
    list_high_keys = []
    list_high_values = []
    for k, v in dict_high.items():
        list_high_keys.append(k)
        list_high_values.append(v)
    # 把覆盖深度低的位置存为列表
    list_low_keys = []
    list_low_values = []
    for k, v in dict_low.items():
        list_low_keys.append(k)
        list_low_values.append(v)
    # 打印出前几行
    #n = 0
    # while n < 10:
        #print(list_high_keys[n], dict_high[list_high_keys[n]])
        #print(list_high_keys[n], list_high_values[n])
        #n += 1

# 折线图
x1 = list_high_keys  # 点1的横坐标
x2 = list_low_keys  # 点2的横坐标
x3 = list_all_keys
k1 = list_high_values  # 线1的纵坐标
k2 = list_low_values  # 线2的纵坐标
k3 = list_all_values
plt.plot(x3, k3, 's-', color='r', label="all coverage")  # s-:方形
# plt.plot(x1, k1, 's-', color='b', label="high coverage")  # s-:方形
plt.plot(x2, k2, 'o-', color='g', label="low coverage")  # o-:圆形

plt.xlabel("pos")  # 横坐标名字
plt.ylabel("depth")  # 纵坐标名字
plt.legend(loc="best")  # 图例
plt.show()
