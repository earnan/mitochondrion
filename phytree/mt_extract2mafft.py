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
    add_help=False, usage='\npython3   mt_from_gbk_get_cds.py\n每个物种都生成cds及完整序列2个文件')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[dir]', help='输入gbk所在目录', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\gbk', required=False)
optional.add_argument('-o', '--output',
                      metavar='[dir]', help='输出的路径,每个物种都生成cds及完整序列2个文件', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\out', required=False)
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


###############################################################
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################
