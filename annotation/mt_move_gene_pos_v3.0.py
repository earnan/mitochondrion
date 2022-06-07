#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_move_gene_pos_v3.0.py
#         Author:   yujie
#    Description:   mt_move_gene_pos_v3.0.py
#        Version:   3.0
#           Time:   2022/06/07 09:38:21
#  Last Modified:   2022/06/07 09:38:21
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   线粒体平移基因修改位置V3.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[file]', type=str, help='要修改的文件,默认', default='gene.annotation.info', required=False)
optional.add_argument(
    '-o', '--output', metavar='[file]', type=str, help='输出文件,默认', default='final_gene.annotation.info', required=False)
optional.add_argument(
    '-ln', '--line_number', metavar='[int]', type=int, help='从第几行开始分开操作,默认不分开', default=1,  required=False)
optional.add_argument(
    '-n1', '--number1', metavar='[int]', type=int, help='平移距离1',  required=False)
optional.add_argument(
    '-n2', '--number2', metavar='[int]', type=int, help='平移距离2',  required=False)
optional.add_argument(
    '-h1', '--info', metavar='[完整的帮助信息]', type=bool, help="使用时'-h1 1'即可", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='帮助信息')
args = parser.parse_args()
if args.info:
    print("\n以下是帮助信息: ")
    print("     适用于内含子外显子")
    print("     #20220210第三版考虑了分段操作,仍需进一步排序\n      分段操作即修改环状的起点\n      此外第三版都改为了子函数")
    print("     增加了排序函数")
    print("\n")


# 替换多个指定位置上的字符(原字符串,位置列表,新的字符列表),使用时位置列表按自然语言习惯

def multi_replace(origin_string, pos_list, new_str_list):  # 把 1-10:- 变为 1-10:minus
    origin_string_list = []
    for s in origin_string:
        origin_string_list.append(s)
    new_string_list = origin_string_list
    # enumerate(sequence, [start=0]) start定义角标从哪开始,默认为0,同索引计数一致
    for i, pos_info in enumerate(pos_list):
        new_string_list[pos_info-1] = new_str_list[i]
    new_string = ''.join(new_string_list)
    return new_string


def edit_pos(pos_info, n):  # 1-10:+ 修改为 101-110:+
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


def get_new_line(line, n):
    line_content = line.split()
    pos_info = line_content[1]  # 例如67574-67687:-;135718-135975:+
    s_content = pos_info.split(';')  # 分号隔开
    if len(s_content) == 1:
        new_pos_info = edit_pos(pos_info, n)
    elif len(s_content) == 2:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        new_pos_info = new_pos_info1+';'+new_pos_info2
    elif len(s_content) == 3:
        new_pos_info1 = edit_pos(s_content[0], n)
        new_pos_info2 = edit_pos(s_content[1], n)
        new_pos_info3 = edit_pos(s_content[2], n)
        new_pos_info = new_pos_info1+';'+new_pos_info2+';'+new_pos_info3
    new_line = line.replace(pos_info, new_pos_info)
    print(new_line.strip('\n'))
    return new_line


# 平移的主函数
fi = open(args.input, 'r')
tmp = open('tmp', 'w')
ln = args.line_number
n1 = args.number1
n2 = args.number2
for i in range(0, ln-1):  # 前(ln-1)行平移n1距离,从第ln行开始平移n2距离
    line = fi.readline()
    new_line = get_new_line(line, n1)
    tmp.write(new_line)
for line in fi:
    new_line = get_new_line(line, n2)
    tmp.write(new_line)
print('\n')
fi.close()
tmp.close()


# 排序的子函数
def pos_sort(input_file, output_file):
    fi = open(input_file, 'r')
    fo = open(output_file, 'w')
    dic = {}
    for line in fi:
        id = int(line.split()[1].split('-')[0])
        dic[id] = line
    # print(dic)
    for i in sorted(dic):
        print(dic[i].strip('\n'))
        fo.write(dic[i])
    return 0


pos_sort('tmp', args.output)
