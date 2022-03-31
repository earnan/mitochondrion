#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_from_gbk_get_cds.py
# Original Author:
#    Description:   mt_from_gbk_get_cds.py
#        Version:   1.0
#           Time:   2022/03/09 15:21:51
#  Last Modified:   2022/03/25 14:05:51
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
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


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


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


def merge_sequence(ele, complete_seq):  # 合并获取到的序列,用于函数(获取cds)来调用
    cds_seq = ""
    tmp_list = []  # 位置列表
    for ele1 in ele.location.parts:
        if ele1.strand == (-1):
            # print('minus')
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际起点,从end中取不用+1
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际终点,从start取+1
            cds_seq += ir(complete_seq[ele1.start:ele1.end])
        elif ele1.strand == (1):
            # print('plus')
            tmp_list.append(str(int(re.findall(
                r'\d+', str(ele1.start))[0])+1))  # 实际起点,要+1
            tmp_list.append(re.findall(
                r'\d+', str(ele1.end))[0])  # 实际终点,不用+1
            # 切片没问题,索引从start到end-1,也就是对应start+1到end的序列
            cds_seq += complete_seq[ele1.start:ele1.end]
    # print(tmp_list)
    return tmp_list, cds_seq


def get_complete_note(seq_record):  # 获取整个完整基因组ID
    seq_id = ''
    # if seq_record.description.find('chloroplast'):#有bug
    if seq_record.description.split(',')[-2].split()[-1] == 'chloroplast':
        seq_id = seq_record.description.split(
            'chloroplast')[0].replace(' ', '_').rstrip('_')
        name = seq_record.name
        if seq_id == name:
            seq_id = seq_id
        elif seq_id != name:
            seq_id = seq_id+'_'+name
        complete_note = ">" + seq_id + "\n"  # chloroplast--叶绿体
    elif seq_record.description.split(',')[-2].split()[-1] == 'mitochondrion':
        seq_id = seq_record.description.split(
            'mitochondrion')[0].replace(' ', '_').rstrip('_')
        name = seq_record.name
        if seq_id == name:
            seq_id = seq_id
        elif seq_id != name:
            seq_id = seq_id+'_'+name
        complete_note = ">" + seq_id + "\n"  # mitochondrion--线粒体
    else:
        print('WARNING')
        complete_note = ">" + (seq_record.description.split('chloroplast')
                               [0]).replace(' ', '_').rstrip('_') + "\n"
    return complete_note, seq_id


def get_cds_note(ele, complete_seq, seq_id):  # 获取cds的id及序列
    if len(ele.location.parts) == 3:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + tmp_list[3]+';' + \
            tmp_list[4]+".." + tmp_list[5]+"]" + " [gene=" + \
            ele.qualifiers['gene'][0] + "]" + "\n"  # '>'后的格式和已有脚本兼容
    elif len(ele.location.parts) == 2:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + \
            tmp_list[3]+"]" + " [gene=" + ele.qualifiers['gene'][0] + \
            "]" + "\n"               # '>'后的格式和已有脚本兼容
    elif len(ele.location.parts) == 1:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+"]" + \
            " [gene=" + ele.qualifiers['gene'][0] + "]" + \
            "\n"    # '>'后的格式和已有脚本兼容
    return cds_note, cds_seq


def get_cds(gbk_file, flag, dict_gene_len):  # 解析gbk文件获取cds
    """完整基因组"""
    seq_record = SeqIO.read(gbk_file, "genbank")
    complete_seq = str(seq_record.seq)
    complete_note, seq_id = get_complete_note(seq_record)
    complete_fasta = format_fasta(complete_note, complete_seq, 70)  # 70换行本例不采用
    """cds序列"""
    count = 0  # 对cds数量计数
    cds_fasta = ""
    list_gene_name = []  # 统计cds
    for ele in seq_record.features:
        if ele.type == "CDS":
            # print(ele)
            # print(ele.qualifiers)
            #print(3*(len(ele.qualifiers['translation'][0])+1), len(cds_seq))
            count += 1
            cds_note, cds_seq = get_cds_note(ele, complete_seq, seq_id)
            cds_fasta += format_fasta(cds_note, cds_seq, 70)

            gene_name = ele.qualifiers['gene'][0]
            gene_name = gene_name_standardization(gene_name)
            list_gene_name.append(gene_name)
            dict_gene_len[gene_name].append(
                3*(len(ele.qualifiers['translation'][0])+1))

            if (flag):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
    file_name = os.path.basename(gbk_file)
    s = '{0}有{1}个CDS'.format(file_name, count)
    print(s)
    print(list_gene_name)
    return cds_fasta, complete_fasta, count, file_name, list_gene_name, s, dict_gene_len, seq_id


def gene_name_standardization(gene_name):  # 格式化基因名字,可重复使用
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            if all_gene_list_lower[i] == gene_name:
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 13:
            print(gene_name)
            print('WARNING!Please check!')
    return gene_name


if __name__ == '__main__':
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    dict_gene_len = {}  # 统计每个基因在不同物种中的长度,取平均
    for i in all_gene_list_upper:
        dict_gene_len[i] = []
    dict_file_cds_count = {}  # 每个文件中cds计数
    dict_missing_gene = {}  # 每个文件中缺失的基因统计
    with open((args.output+os.sep+'log'), 'w') as f_log:
        f_log.write('gene{0}atp6{0}atp8{0}cob{0}cox1{0}cox2{0}cox3{0}nad1{0}nad2{0}nad3{0}nad4{0}nad4L{0}nad5{0}nad6\n'.format(
            '\t'))
        f_log.write('gene{0}ATP6{0}ATP8{0}CYTB{0}COX1{0}COX2{0}COX3{0}ND1{0}ND2{0}ND3{0}ND4{0}ND4L{0}ND5{0}ND6\n'.format(
            '\t'))

    file_list = os.listdir(args.input)
    file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字
    # print(file_list)
    for file in file_list:
        ingbk_path = os.path.join(args.input, file)
        cds_fasta, complete_fasta, count, file_name,  list_gene_name, s, dict_gene_len, seq_id = get_cds(
            ingbk_path, False, dict_gene_len)  # , all_gene_list_upper, all_gene_list_lower)
        # dict_file_cds_count[file_name] = count  # 每个文件中cds计数
        dict_file_cds_count[seq_id] = count  # 每个文件中cds计数

        with open((args.output+os.sep+file_name.rstrip('.gbk')+'_complete.fasta'), 'w') as f_complete, \
                open((args.output+os.sep+file_name.rstrip('.gbk')+'_cds.fasta'), 'w') as f_cds, \
                open((args.output+os.sep+'log'), 'a+') as f_log:
            f_cds.write(cds_fasta)
            f_complete.write(complete_fasta)
            f_log.write(s+'\n')
            f_log.write('>'+file.rstrip('.gbk')+'\t')
            list_missing_gene = []
            for i in range(len(all_gene_list_upper)):
                if all_gene_list_upper[i] in list_gene_name \
                        or all_gene_list_upper[i].lower() in list_gene_name \
                        or all_gene_list_lower[i] in list_gene_name \
                        or all_gene_list_lower[i].upper() in list_gene_name:
                    f_log.write(all_gene_list_upper[i]+'\t')
                else:
                    f_log.write('NULL'+'\t')
                    list_missing_gene.append(all_gene_list_upper[i])  # 缺失的基因
                    # print(sum(dict_gene_len[all_gene_list_upper[i]])/len(dict_gene_len[all_gene_list_upper[i]]))
            print(list_missing_gene)
            f_log.write('\n')
            [f_log.write(tmp+'\t') for tmp in list_missing_gene]
            f_log.write('\n')
        dict_missing_gene['>'+seq_id] = list_missing_gene
    with open((args.output+os.sep+'log'), 'a+') as f_log:
        f_log.write(str(dict_missing_gene))

    total_ref_gene = 0  # 除了物种1外,其他物种所有cds总数
    for i in dict_file_cds_count.keys():
        total_ref_gene += dict_file_cds_count[i]
    print(dict_gene_len)
    print(dict_missing_gene)
    print(2*(total_ref_gene-13))
    print('\n')
#######################
    # 用gap构造没有的基因
    if args.check:
        for i in dict_missing_gene.keys():
            cds_fasta = ''
            for j in dict_missing_gene[i]:
                # print(j)
                ave = round(sum(dict_gene_len[j]) /
                            len(dict_gene_len[j]))  # 该基因平均长度
                cds_note = (i+' [0..0]'+' [gene={}]').format(j)
                cds_seq = ave*'-'
                cds_fasta += format_fasta(cds_note, cds_seq, 70)
            print(cds_fasta)
            file_name = (i.split('_')[-2]+'_' +
                         i.split('_')[-1]+'.1').lstrip('>')
            with open(args.output+os.sep+file_name+'_cds.fasta', 'a+') as f_cds:
                f_cds.write(cds_fasta)

###############################################################
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
print('Done')
###############################################################