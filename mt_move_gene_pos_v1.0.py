#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename: 编程处理模板.py
#         Author: yuj@genepioneer.cn
#    Description: sample
#  Last Modified: 2021-xx-xx 16:29:29
#
# Copyright (C) 2021xxxx genepioneer Corporation
##########################################################
import argparse

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   线粒体平移基因修改位置(只适用不分段基因)')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[file]', type=str, help='要修改的文件,默认', default='gene.annotation.info', required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='输出文件,默认', default='final_gene.annotation.info', required=False)
optional.add_argument(
    '-ln', '--line_number', metavar='[int]', type=int, help='从第ln行开始,默认1', default=1,  required=False)
optional.add_argument(
    '-n', '--number', metavar='[int]', type=int, help='平移距离',  required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()


# 替换多个指定位置上的字符(原字符串,位置列表,新的字符列表),使用时位置列表按自然语言习惯
def multi_replace(origin_string, pos_list, new_str_list):
    origin_string_list = []
    for s in origin_string:
        origin_string_list.append(s)
    new_string_list = origin_string_list
    # enumerate(sequence, [start=0]) start定义角标从哪开始,默认为0,同索引计数一致
    for i, pos_info in enumerate(pos_list):
        new_string_list[pos_info-1] = new_str_list[i]
    new_string = ''.join(new_string_list)
    return new_string


fi = open(args.input, 'r')
fo = open(args.output, 'w')
ln = args.line_number
n = args.number
for i in range(0, ln-1):  # 从第ln行开始也就是跳过ln-1行
    line = fi.readline()
    fo.write(line)
for line in fi:
    line_content = line.split()
    pos_info = line_content[1]  # 例如1-10:+

    if pos_info.endswith('-'):
        pos_info = multi_replace(pos_info, [0], ['minus'])  # 索引-1对应自然的0

    pos_content = pos_info.split('-')
    new_start = int(pos_content[0])+n  # 例如1

    end_pos_all = pos_content[1]  # 例如10:+
    end_pos_content = end_pos_all.split(':')
    new_end = int(end_pos_content[0])+n  # 例如10
    new_pos_info = (str(new_start)+'-'+str(new_end)+':' +
                    end_pos_content[1]).replace('minus', '-')

    new_line = line.replace(pos_info, new_pos_info)
    print(new_line)
    fo.write(new_line)

fi.close()
fo.close()
# 测试
