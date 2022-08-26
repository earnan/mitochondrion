#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_from_gbk_get_cds.py
#         Author:   yujie
#    Description:   mt_from_gbk_get_cds.py
#        Version:   1.0
#           Time:   2022/03/09 15:21:51
#  Last Modified:   2022/08/19 11:59:10
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
#import pretty_errors
import re
import sys
import time
#import copy
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   mt_from_gbk_get_cds_V2.5.py\n\
#         Author:   yujie\n\
#    Description:   mt_from_gbk_get_cds_V2.5.py\n\
#        Version:   2.0\n\
#           Time:   2022/03/09 15:21:51\n\
#  Last Modified:   2022/08/19 11:59:10\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   mt_from_gbk_get_cds_V1.0.py\n\
功能：每个物种都生成cds及完整序列2个文件\n\
1.常规使用\n\
1.1 -i [gbk dir] -o [out dir] \n\
2.其他使用\n\
2.1 -c 用GAP构造缺失的cds\n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\mitochondrion\phytree\mt_from_gbk_get_cds_V1.0.py\n\
Path: /share/nas1/yuj/script/mitochondrion/phytree/mt_from_gbk_get_cds_V1.0.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[indir]', help='gbk dir', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\gbk',  required=False)
optional.add_argument(
    '-o', '--output', metavar='[outdir]', help='输出的路径', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\cds', required=False)
optional.add_argument('-c', '--check', help='默认否,使用时-c',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

# ############################################################################################


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note.strip() + "\n" + format_seq + "\n"


def ir(s):  # 反向互补
    re = s[::-1]  # 字符串反向
    c = ""
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
    try:
        seq_id = ''
        # 20220819 NC_044756.1.gbk voucher Liu HM/CP02 chloroplast  有特殊符号，需要处理
        seq_record.description = seq_record.description.replace('/', '_')
        # if seq_record.description.find('chloroplast'):#有bug,用str格式化后就没问题了
        # 20220627 if str(seq_record.description).find('chloroplast') -1也成立,判断时一定要以False True为准
        # or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' or seq_record.description.split(',')[-2].split()[-1] == 'plastid':
        if str(seq_record.description).find('chloroplast') > 0:
            seq_id = seq_record.description.split(
                'chloroplast')[0].replace(' ', '_').rstrip('_')
            name = seq_record.name
            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"
        # or seq_record.description.split(',')[-2].split()[-1] == 'mitochondrion':
        # 20220811 NC_031548.gbk 描述部分 DEFINITION  Oxynoemacheilus angorae mitochondrial DNA, complete genome.
        elif str(seq_record.description).find('mitochondrion') > 0 or str(seq_record.description).find('mitochondrial') > 0:
            seq_id = seq_record.description.split(
                'mitochondrion')[0].replace(' ', '_').rstrip('_')  # 物种或样品名

            if seq_id.startswith('UNVERIFIED:_'):  # 去掉 UNVERIFIED:_
                seq_id = seq_id.lstrip('UNVERIFIED:_')
            # 去掉Cerion_watlingense_voucher_USNM:1514170_MN904501 中 冒号后的内容
            if len(seq_id.split(':')) > 1:
                seq_id = seq_id.split(':')[0]

            name = seq_record.name  # 要么是登录号  要么是样本

            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"
        else:
            print('Genome Type WARNING! {}!'.format(
                seq_record.description.split(', ')[-2].split()[-1]))
            complete_note = ">" + \
                (seq_record.description.split('chloroplast')
                 [0]).replace(' ', '_').rstrip('_') + "\n"
    except:  # 如果遇到任何出错
        print('try/except')
        complete_note = ''
        #gbk_type = input('genome type(1:chloroplast;2:mitochondrion): ')
        gbk_type = 2
        if gbk_type == 2:
            seq_id = seq_record.description.split(
                'mitochondrion')[0].replace(' ', '_').rstrip('_')  # 物种或样品名

            if seq_id.startswith('UNVERIFIED:_'):  # 去掉 UNVERIFIED:_
                seq_id = seq_id.lstrip('UNVERIFIED:_')
            # 去掉Cerion_watlingense_voucher_USNM:1514170_MN904501 中 冒号后的内容
            if len(seq_id.split(':')) > 1:
                seq_id = seq_id.split(':')[0]

            name = seq_record.name  # 要么是登录号  要么是样本

            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"

    return complete_note, seq_id


# 仅在 get_cds_note(ele, complete_seq, seq_id, tmp_gene_name)中使用，以应对"gbk中cds没有/gene标签，但是有/product标签"的情况
def gene_name_standardization_1(gene_name):  # 格式化基因名字,可重复使用
    #name_flag = 0
    all_gene_dict = {'ATP synthase F0 subunit 6': 'ATP6', 'ATP synthase F0 subunit 8': 'ATP8', 'cytochrome b': 'CYTB',
                     'cytochrome c oxidase subunit I': 'COX1', 'cytochrome c oxidase subunit II': 'COX2',
                     'cytochrome c oxidase subunit III': 'COX3', 'NADH dehydrogenase subunit 1': 'ND1',
                     'NADH dehydrogenase subunit 2': 'ND2', 'NADH dehydrogenase subunit 3': 'ND3',
                     'NADH dehydrogenase subunit 4': 'ND4', 'NADH dehydrogenase subunit 4L': 'ND4L',
                     'NADH dehydrogenase subunit 5': 'ND5', 'NADH dehydrogenase subunit 6': 'ND6'}
    if gene_name in all_gene_dict.keys():
        gene_name = all_gene_dict[gene_name]
    '''
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    gene_name = gene_name.replace('III', '3').replace(
        'II', '2').replace('I', '1')  # 20220825 gbk文件里 基因名是罗马数字，因此需要先处理
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            # 20220624  JN619347.gbk 有一个nad4L 大小写混合 就离谱
            if all_gene_list_lower[i] == gene_name.lower():
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 13:
            print(gene_name, 'WARNING!Please check!')
            name_flag = 1
    '''
    return gene_name  # , name_flag


def get_cds_note(ele, complete_seq, seq_id, tmp_gene_name):  # 获取cds的id及序列

    # 20220825 NC_034226.gbk cds没有/gene标签，但是有/product标签
    if 'gene' not in ele.qualifiers.keys():
        if 'product' in ele.qualifiers.keys():
            tmp_gene_name = ele.qualifiers['product'][0]
            tmp_gene_name = gene_name_standardization_1(tmp_gene_name)
        else:
            tmp_gene_name = input(
                "\nPrevious: [{0}]. Current: {1}.\nPlease input current gene name:".format(tmp_gene_name, ele.location.parts))  # 返回上一个基因,好从其他参考找这个没名字的
        '''
        try:
            #tmp_gene_name = ele.qualifiers['note'][0]
            print(ele.qualifiers)
        except:
            tmp_gene_name = input(
                "\nPrevious: [{0}]. Current: {1}.\nPlease input current gene name:".format(tmp_gene_name, ele.location.parts))
        '''
    elif 'gene' in ele.qualifiers.keys():
        tmp_gene_name = ele.qualifiers['gene'][0]

    if len(ele.location.parts) == 3:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + tmp_list[3]+';' + \
            tmp_list[4]+".." + tmp_list[5]+"]" + " [gene=" + \
            tmp_gene_name + "]" + "\n"  # '>'后的格式和已有脚本兼容

    elif len(ele.location.parts) == 2:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+';' + tmp_list[2]+".." + \
            tmp_list[3]+"]" + " [gene=" + tmp_gene_name + \
            "]" + "\n"               # '>'后的格式和已有脚本兼容

    elif len(ele.location.parts) == 1:
        tmp_list, cds_seq = merge_sequence(ele, complete_seq)
        cds_note = ">" + seq_id + " [" + tmp_list[0]+".." + tmp_list[1]+"]" + \
            " [gene=" + tmp_gene_name + "]" + \
            "\n"    # '>'后的格式和已有脚本兼容

    return cds_note, cds_seq, tmp_gene_name


def gene_name_standardization_2(gene_name):  # 格式化基因名字,可重复使用
    name_flag = 0
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    gene_name = gene_name.replace('III', '3').replace(
        'II', '2').replace('I', '1')  # 20220825 gbk文件里 基因名是罗马数字，因此需要先处理
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            # 20220624  JN619347.gbk 有一个nad4L 大小写混合 就离谱
            if all_gene_list_lower[i] == gene_name.lower():
                gene_name = all_gene_list_upper[i]
                break
            else:
                i += 1
        if i >= 13:
            print(gene_name, 'WARNING!Please check!')
            name_flag = 1
    return gene_name, name_flag


def get_cds(gbk_file, flag, dict_gene_len, file_no):  # 解析gbk文件获取cds
    """完整基因组"""
    seq_record = SeqIO.read(gbk_file, "genbank")
    complete_seq = str(seq_record.seq)

    complete_note, seq_id = get_complete_note(seq_record)
    complete_fasta = format_fasta(complete_note, complete_seq, 70)  # 70换行本例不采用

    """cds序列"""
    count = 0  # 对cds数量计数
    cds_fasta = ""
    list_gene_name = []  # 统计cds
    tmp_gene_name = ''  # 上一个基因名字,为子函数get_cds_note()准备的
    for ele in seq_record.features:
        if ele.type == "CDS":
            count += 1
            cds_note, cds_seq, tmp_gene_name = get_cds_note(
                ele, complete_seq, seq_id, tmp_gene_name)
            # list_gene_name.append(tmp_gene_name)  # 本次的基因名字 复用,线粒体的话,在下一部分存入列表
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            gene_name = tmp_gene_name
            gene_name, name_flag = gene_name_standardization_2(gene_name)
            if name_flag == 1:
                print(gbk_file, 'WARNING!Please check!')
            list_gene_name.append(gene_name)  # 存入列表
            dict_gene_len[gene_name].append(
                3*(len(ele.qualifiers['translation'][0])+1))  # cds序列长度
            if (flag):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
    file_name = os.path.basename(gbk_file)
    s = '{2}: {0}有{1}个CDS'.format(file_name, count, file_no)

    if count == 0:
        # --------There may be no comments--------'.format(
        s = '{2}: {0}有{1}个CDS'.format(file_name, count, file_no)
        print(s.ljust(50), '----------There may be no comments----------')
    elif count != 0:
        # 20220811 输出左对齐  str.ljust(50)  达到50个字符宽度
        print(s.ljust(50), '+', list_gene_name)

    return cds_fasta, complete_fasta, count, file_name, list_gene_name, s, dict_gene_len, seq_id


def create_gene_by_gap(dict_missing_gene, dict_gene_len, cds_file_path):  # 用gap构造没有的基因
    for i in dict_missing_gene.keys():
        cds_fasta = ''
        for j in dict_missing_gene[i]:
            ave = round(sum(dict_gene_len[j]) /
                        len(dict_gene_len[j]))  # 该基因平均长度
            cds_note = (i+' [0..0]'+' [gene={}]').format(j)
            cds_seq = ave*'-'
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
        print(cds_fasta)
        file_name = (i.split('_')[-2]+'_' +
                     i.split('_')[-1]+'.1').lstrip('>')
        with open(cds_file_path, 'ab+') as f_cds:  # 读写打开一个二进制文件，允许读，或在文件末追加数据
            f_cds.write(cds_fasta.encode())


if __name__ == '__main__':
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################
    print('\n')

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    """写入统计文件"""
    with open((args.output+os.sep+'log'), 'w') as f_log:
        f_log.write(
            'gene{0}atp6{0}atp8{0}cob{0}cox1{0}cox2{0}cox3{0}nad1{0}nad2{0}nad3{0}nad4{0}nad4L{0}nad5{0}nad6\n'.format('\t'))  # 小写
        f_log.write(
            'gene{0}ATP6{0}ATP8{0}CYTB{0}COX1{0}COX2{0}COX3{0}ND1{0}ND2{0}ND3{0}ND4{0}ND4L{0}ND5{0}ND6\n'.format('\t'))  # 大写
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    """统计初始化"""
    dict_missing_gene = {}  # 每个文件中缺失的基因统计,总 字典
    dict_gene_len = {}  # 统计每个基因在不同物种中的长度,取平均
    for i in all_gene_list_upper:
        dict_gene_len[i] = []

    """初始化"""
    dict_file_cds_count = {}  # 每个文件中cds计数
    file_list = [i for i in os.listdir(
        args.input) if i.endswith('gbk')]  # 20220811 检查文件后缀
    file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字

    """主程序"""
    file_no = 0
    for file in file_list:
        file_no += 1
        ingbk_path = os.path.join(args.input, file)
        cds_fasta, complete_fasta, count, file_name,  list_gene_name, s, dict_gene_len, seq_id = get_cds(
            ingbk_path, False, dict_gene_len, file_no)
        dict_file_cds_count[seq_id] = count  # 每个文件中cds计数
        """写入文件"""
        with open((args.output+os.sep+seq_id+'.fasta'), 'wb') as f_complete, \
                open((args.output+os.sep+'cds_'+file_name.rstrip('.gbk')+'.fasta'), 'wb') as f_cds, \
                open((args.output+os.sep+'log'), 'a+') as f_log:
            f_cds.write(cds_fasta.encode())
            f_complete.write(complete_fasta.encode())

            """以下为统计部分"""
            f_log.write(s+'\n')
            f_log.write('>'+file.rstrip('.gbk')+'\t')
            list_missing_gene = []  # 每个文件中缺失的基因统计,单列表
            for i in range(len(all_gene_list_upper)):
                if all_gene_list_upper[i] in list_gene_name \
                        or all_gene_list_upper[i].lower() in list_gene_name \
                        or all_gene_list_lower[i] in list_gene_name \
                        or all_gene_list_lower[i].upper() in list_gene_name:
                    f_log.write(all_gene_list_upper[i]+'\t')
                else:
                    f_log.write('NULL'+'\t')
                    list_missing_gene.append(all_gene_list_upper[i])  # 缺失的基因
            if 1 <= len(list_missing_gene) < 13:
                print(''.ljust(55), '-', list_missing_gene)

            f_log.write('\n')
            [f_log.write(tmp+'\t') for tmp in list_missing_gene]
            f_log.write('\n')
        dict_missing_gene['>'+seq_id] = list_missing_gene  # 放入字典当中
    with open((args.output+os.sep+'log'), 'a+') as f_log:
        f_log.write(str(dict_missing_gene))

    """除了物种1外,其他物种所有cds总数"""
    total_ref_gene = 0
    for i in dict_file_cds_count.keys():  # 键为seq_id,值为个数
        total_ref_gene += dict_file_cds_count[i]
    print('\n')
    print(dict_gene_len)  # 键为每个基因,值为列表,列表为每个基因在不同物种中的长度
    print(dict_missing_gene)  # 键为>seq_id,值为列表,列表为每个物种确实的基因
    print(2*(total_ref_gene-13))

    """gap构造基因"""
    if args.check:
        cds_file_path = args.output+os.sep+file_name+'_cds.fasta'
        create_gene_by_gap(dict_missing_gene, dict_gene_len, cds_file_path)

    print('\n')
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
