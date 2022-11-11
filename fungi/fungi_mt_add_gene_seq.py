#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   fungi_mt_add_gene_seq.py
#         Author:   yujie
#    Description:   fungi_mt_add_gene_seq.py
#        Version:   2.0
#           Time:   2022/11/09 14:19:17
#  Last Modified:   2022/11/09 14:19:17
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from humre import *  # 正则
# from icecream import ic  # 打印
import argparse  # 命令行
# import linecache  # 大文件行读取
import os  # 目录路径
# import pretty_errors  # 错误提示
import re  # 正则
import sys
# import time
# import copy  # 深度拷贝
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   fungi_mt_add_gene_seq.py\n\
#         Author:   yujie\n\
#    Description:   fungi_mt_add_gene_seq.py\n\
#        Version:   2.0\n\
#           Time:   2022/11/09 16:35:19\n\
#  Last Modified:   2022/11/09 15:52:22\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   GNU General Public License v3.0\n\
#\n\
##########################################################\n\
\n\
\npython3   fungi_mt_add_gene_seq.py\n\
\n\
1.常规使用\n\
1.1查看密码子 -n -i -p \n\
1.2序列被误判为trna,强制翻译使用 -f \n\
\n\
2.递归查找与存储\n\
2.1起始子自动查找,-m 最大查找次数 -df 向后查找（默认向前查找）\n\
2.2终止子自动查找,-m 最大查找次数\n\
\n\
3.存储为文件\n\
3.1存储序列,-sn 基因名\n\
3.2存储蛋白,-sp 基因名\n\
\n\
Path: E:\OneDrive\jshy信息部\Script\mitochondrion\fungi\mt_add_gene_seq.py\n\
Path: /share/nas1/yuj/script/mitochondrion/fungi/mt_add_gene_seq.py\n\
Version: V2.0'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infasta', metavar='[infasta]', help='输入fasta文件', type=str, default='E:\\OneDrive\\jshy信息部\\Script\\mitochondrion\\fungi\\Ustilago_esculenta_FULLMT.fsa', required=False)
optional.add_argument(
    '-p', '--pos_str', metavar='[pos_str]', help="输入位置,形如'124353-124892:-;126001-126552:-'", type=str, default='6130-6378:-', required=False)

optional.add_argument(
    '-n', '--codonnumber', metavar='[codon_number]', help='密码子表,默认4', type=int, default=4, required=False)
optional.add_argument(
    '-m', '--maxnumber', metavar='[max_number]', help='最大递归查找次数,默认0,假查找', type=int, default=0, required=False)
optional.add_argument('-df', '--direction_flag',
                      help='起始子查找方向,默认true向前(序列变长),向后则-df', action='store_false', required=False)
optional.add_argument('-trans', '--trans_flag',
                      help='翻译?默认是,不运行则-trans', action='store_false', required=False)
optional.add_argument('-f', '--force_flag',
                      help='强制翻译?默认否,运行则-f', action='store_true', required=False)
optional.add_argument('-sn', '--nuc_file_name',
                      metavar='[store 2 dna]', help='默认否,值为NULL,存储则输入gene名', type=str,  default='NULL', required=False)
optional.add_argument('-sp', '--pro_file_name',
                      metavar='[store 2 protein]', help='默认否,值为NULL,存储则输入蛋白名', type=str,  default='NULL', required=False)
optional.add_argument('-info', help='更新日志,使用时-info',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20221110  第一版，待解决递归退出问题')
    print('\n')
    sys.exit(0)

##########################################################################################


def read_file(infasta):  # 读取文件
    with open(infasta, 'r') as f:
        full_seq = ''
        for line in f:
            if not line.startswith('>'):
                full_seq += line.strip('\n')
    return full_seq


def ir(s):  # 反向互补
    re = s[::-1]
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


def format_pos(pos_info):  # 读取输入的位置为位置列表 或 反过来转化
    pos_list_info = []
    pos_str_info = ''
    if type(pos_info) == type(''):
        content = pos_info.split(';')
        for ele in content:
            if ele.split(':')[-1] == '-':
                tmp = ele.split(':')[0]+':'+'-1'
                pos_list_info.append(tmp)
            elif ele.split(':')[-1] == '+':
                tmp = ele.split(':')[0]+':'+'1'
                pos_list_info.append(tmp)
        pos_info = pos_list_info
    elif type(pos_info) == type([]):
        for ele in pos_list_info:
            if ele.split(':')[-1] == '-1':
                tmp = ele.split(':')[0]+':'+'-'
                pos_str_info += (tmp+';')
            elif ele.split(':')[-1] == '1':
                tmp = ele.split(':')[0]+':'+'+'
                pos_str_info += (tmp+';')
        pos_info = pos_str_info.rstrip(';')
    return pos_info

#######################################################################################################################


def merge_sequence(pos_info, full_seq):  # 合并获取到的序列 并排序位置 排序：把位置列表按序列方向排
    if type(pos_info) == type([]):
        pos_list = pos_info
    """
    20220728 判断是否是trna,返回一个flag
    """
    flag_gene_type = 'NULL'
    len_trna_gene = 0
    if type(pos_info) == type([]) and len(pos_info) == 1:  # 位置列表只有一段，直接计算长度即可

        start = re.findall(r'[0-9]+', pos_info[0].split(':')[0])[0]
        end = re.findall(r'[0-9]+', pos_info[0].split(':')[0])[-1]
        len_trna_gene = abs(int(end)-int(start))+1
        if 55 <= len_trna_gene <= 100:
            flag_gene_type = 'trna'

    """
    20221101 参考线粒体  解决跨首尾基因
    """
    if type(pos_info) == type([]) and len(pos_info) == 1 and pos_info[0].split(':')[-1] == '1' and \
        int(re.findall(r'[0-9]+', pos_info[0].split(':')[0])[0]) > \
            int(re.findall(r'[0-9]+', pos_info[0].split(':')[0])[-1]):  # 14323-1527:1
        pos1 = '{0}-{1}:1'.format(re.findall(r'[0-9]+',
                                  pos_info[0].split(':')[0])[0], len(full_seq))
        pos2 = '1-{}:1'.format(re.findall(r'[0-9]+',
                               pos_info[0].split(':')[0])[-1])
        pos_list = [pos1, pos2]

    cds_seq = ""
    if int(pos_list[0].split(':')[-1]) == -1:  # 负链的话
        sorted_pos_list = pos_list[::-1]
    else:
        sorted_pos_list = pos_list  # 正向的话，已经是排序的

    for ele in sorted_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        start_index = start-1
        end_index = end
        if strand == (-1):
            # seq[start_index:end_index] 角标从start_index到end_index
            # 取的是索引start-1一直到end  取的是start一直到end的碱基
            cds_seq += ir(full_seq[start_index:end_index])
        elif strand == (1):
            cds_seq += full_seq[start_index:end_index]
    return cds_seq, sorted_pos_list, flag_gene_type, len_trna_gene
#################################################################################################################

# 20220722 新增子函数


def storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq):  # 存储获取到的dna序列或蛋白
    # 20220629   trna 存起来
    if flag_gene_type == 'trna' and (not args.force_flag):
        print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
    # 20221101 和线粒体一样 精简
    current_abs_path = os.getcwd()
    if nuc_file_name != 'NULL':
        with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
            f_handle.write(cds_seq+'\n')
    return 0
#######################################################################################################################


def trans2protein_seq(cds_seq, n):  # 翻译成氨基酸,返回是否正确以及第一个终止子在基因序列上的相对位置
    if n == 5:
        start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 5
    elif n == 2:
        start_codon_list = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'AGA',
                          'AGG', 'TA', 'T', 'AG']  # 2
    elif n == 4:
        start_codon_list = ['TTA', 'TTG', 'CTG',
                            'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 4

    seq_check_flag = 0  # seq_check_flag 起始是否正确的标志,默认False   20220610改为数字,0为正确,1为起始x,2为内部错,3为末尾错
    first_stop_codon_index_in_protein_seq = 0
    if len(cds_seq) % 3 == 1:
        print('current len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
    elif len(cds_seq) % 3 == 2:
        print('current len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))

    coding_dna = Seq(cds_seq)
    protein_seq = coding_dna.translate(table=n)
    print('------------------------------------------------------------')
    print(protein_seq)
    if not cds_seq[0:3] in start_codon_list:  # 起始不正确
        print('------------------------------------------------------------#####start is wrong!')
        seq_check_flag = 1
    else:  # 起始正确,又分3种情况
        if protein_seq.count('*') > 1:  # 终止多于1,意味着提前终止
            print(
                '------------------------------------------------------------#####interior is wrong!')
            seq_check_flag = 2
            first_stop_codon_index_in_protein_seq = protein_seq.find('*')
            print('Index of the first stop codon :{}'.format(
                first_stop_codon_index_in_protein_seq))
            print('\n')
        elif protein_seq.count('*') < 1:  # 终止小于1 1.真的未终止2.线粒体终止了,共5种细分情况
            if len(cds_seq) % 3 == 1 and cds_seq[-1] in end_codon_list:
                print(
                    '------------------------------------------------------------ok')
            elif len(cds_seq) % 3 == 2 and cds_seq[-2:] in end_codon_list:
                print(
                    '------------------------------------------------------------ok')
            else:
                print(
                    '------------------------------------------------------------#####end is wrong!')
                seq_check_flag = 3
        elif protein_seq.count('*') == 1:  # 1个终止子,提前终止或者没问题
            if not protein_seq.endswith('*'):
                print(
                    '------------------------------------------------------------#####interior is wrong!')
                seq_check_flag = 2
                first_stop_codon_index_in_protein_seq = protein_seq.find('*')
                print('Index of the first stop codon :{}'.format(
                    first_stop_codon_index_in_protein_seq))
            else:
                seq_check_flag = 0
                print('------------------------------------------------------------ok')
    return seq_check_flag, first_stop_codon_index_in_protein_seq, protein_seq


###################################################################################################################

# 循环查找   *.fas/"1-10:-;20-30:-"/翻译/递归计数/最大递归次数
def loop_look(direction_flag, infasta, pos_str, trans_flag, loop_count, maxnumber, n, nuc_file_name, pro_file_name):
    if n == 5:
        start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 5
    elif n == 2:
        start_codon_list = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'AGA',
                          'AGG', 'TA', 'T', 'AG']  # 2,转录时要加A
    elif n == 4:
        start_codon_list = ['TTA', 'TTG', 'CTG',
                            'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 4
    first_stop_codon_index_in_protein_seq = 0
    modified_pos_str = pos_str  # 最后要传出来，给个初始化
    seq = read_file(infasta)
    pos_list = format_pos(pos_str)
    cds_seq, sorted_pos_list, flag_gene_type, len_trna_type = merge_sequence(
        pos_list, seq)  # sorted_pos_list  把位置当列表再传出来,这个位置信息向下传递
    print('\n'+cds_seq)
    if maxnumber > 0:
        print('current pos:{}'.format(sorted_pos_list))
        print('current pos:{}'.format(pos_str))
    storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq)
    if (trans_flag and (flag_gene_type != 'trna')) or args.force_flag:  # 翻译
        seq_check_flag, first_stop_codon_index_in_protein_seq, protein_seq = trans2protein_seq(
            cds_seq, n)
        if maxnumber > 0:
            if seq_check_flag != 1:
                print('\n[START CONDON]The correct start codon was found after {} searches / Total times: {}'.format(
                    len(loop_count_flag1), len(loop_count_flag)))
                if seq_check_flag != 3:
                    print('\n[STOP CONDON]The correct stop codon was found after {} searches / Total times: {}'.format(
                        len(loop_count_flag3), len(loop_count_flag)))
        current_abs_path = os.getcwd()
        if pro_file_name != 'NULL':
            with open(os.path.join(current_abs_path, pro_file_name+'.protein_seq'), 'w') as f_handle:
                f_handle.write(str(protein_seq)+'\n')
        # #####################################################################################################################################
        # 第一层if else
        # 0为正确,1为起始x,2为内部错,3为末尾错
        if seq_check_flag == 1:  # 起始错,这个优先要满足的条件就不对
            loop_count += 1
            loop_count_flag.append(000)
            loop_count_flag1.append(1)  # 每有一次查找,列表元素个数就+1
            if maxnumber > 0:
                print(
                    '[START CONDON]Start search......Times:{} / Total times:{}'.format(loop_count, len(loop_count_flag)))

            '''
            20221101 叶绿体
            分向前找(flase)  向后找(args.direction_flag=true)
            '''
            if direction_flag == False:
                '''起始向后找 长度减小'''
                # ##############################################################
                # 定义为第二层if else
                # 20220805  如果为假查找，就不进行下一步了
                if cds_seq[3:6] not in start_codon_list and maxnumber > 0:
                    print('old pos:{}'.format(sorted_pos_list))
                    print('old pos:{}'.format(pos_str))
                    # 20220808 以下自动返回位置，也就是开头往后挪6bp
                    if pos_str.split(':')[-1] == '+':
                        modified_pos_str = pos_str.replace(pos_str.split(
                            '-')[0], str(int(pos_str.split('-')[0])+6))
                    elif pos_str.split(':')[-1] == '-':
                        modified_pos_str = pos_str.replace(re.findall(
                            r'\d+', pos_str)[-1], str(int(re.findall(
                                r'\d+', pos_str)[-1])-6))

                    if loop_count <= maxnumber:
                        loop_look(direction_flag, infasta, modified_pos_str, trans_flag,
                                  loop_count, maxnumber, n, nuc_file_name, pro_file_name)
                elif cds_seq[3:6] in start_codon_list:
                    print('\n'+cds_seq)
                    # 20220808 以下自动返回正确位置，也就是开头往后挪3bp
                    if pos_str.split(':')[-1] == '+':
                        modified_pos_str = pos_str.replace(pos_str.split(
                            '-')[0], str(int(pos_str.split('-')[0])+3))
                    elif pos_str.split(':')[-1] == '-':
                        modified_pos_str = pos_str.replace(re.findall(
                            r'\d+', pos_str)[-1], str(int(re.findall(
                                r'\d+', pos_str)[-1])-3))
                    seq_check_flag, first_stop_codon_index_in_protein_seq, protein_seq = trans2protein_seq(
                        cds_seq, n)

                    if maxnumber > 0:
                        if seq_check_flag != 1:
                            print('The correct starting codon was found after {} searches / Total times: {}'.format(
                                len(loop_count_flag1), len(loop_count_flag)))
                    if seq_check_flag == 0:  # 20221020 其他地方还会出错,所以要再次检查
                        print('Correct Position: [{}]'.format(
                            modified_pos_str))
                        if pro_file_name != 'NULL':
                            with open(os.path.join(current_abs_path, pro_file_name+'.protein_seq'), 'w') as f_handle:
                                f_handle.write(str(protein_seq)+'\n')
                    else:  # 还没完全找对
                        if loop_count <= maxnumber:
                            loop_look(direction_flag, infasta, modified_pos_str, trans_flag, loop_count,
                                      maxnumber, n, nuc_file_name, pro_file_name)

            elif direction_flag == True:
                '''起始向前找 长度增加'''
                if maxnumber > 0:
                    print('old pos:{}'.format(sorted_pos_list))
                    print('old pos:{}'.format(pos_str))
                    if pos_str.split(':')[-1] == '+':
                        modified_pos_str = pos_str.replace(pos_str.split(
                            '-')[0], str(int(pos_str.split('-')[0])-3))
                    elif pos_str.split(':')[-1] == '-':
                        modified_pos_str = pos_str.replace(re.findall(
                            r'\d+', pos_str)[-1], str(int(re.findall(
                                r'\d+', pos_str)[-1])+3))

                    if loop_count <= maxnumber:
                        loop_look(direction_flag, infasta, modified_pos_str, trans_flag,
                                  loop_count, maxnumber, n, nuc_file_name, pro_file_name)

        # ################################################################################################################################
        # 第一层if else
        # 0为正确,1为起始x,2为内部错,3为末尾错
        elif seq_check_flag == 3:  # 3代表未终止 序列长度为3种情况,3n 3n+1 3n+2
            loop_count_flag.append(000)
            loop_count_flag3.append(3)
            loop_count = len(loop_count_flag)
            if maxnumber > 0:
                print(
                    '[STOP CONDON]Start search......Times:{} / Total times:{}'.format(len(loop_count_flag3), len(loop_count_flag)))
                print('old pos:{}'.format(sorted_pos_list))
                print('old pos:{}'.format(pos_str))
            if len(cds_seq) % 3 == 0:
                if maxnumber > 0:
                    print(
                        'old len(sequence) is a multiple of three! {}=3n'.format(len(cds_seq)))
                if pos_str.split(':')[-1] == '+':
                    modified_pos_str = pos_str.replace(pos_str.split(
                        '-')[-1].split(':')[0], str(int(pos_str.split('-')[-1].split(':')[0])+3))
                elif pos_str.split(':')[-1] == '-':
                    modified_pos_str = pos_str.replace(re.findall(
                        r'\d+', pos_str)[0], str(int(re.findall(
                            r'\d+', pos_str)[0])-3))
            else:
                if len(cds_seq) % 3 == 1:
                    if maxnumber > 0:
                        print(
                            'old len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
                elif len(cds_seq) % 3 == 2:
                    if maxnumber > 0:
                        print(
                            'old len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))
                if pos_str.split(':')[-1] == '+':
                    modified_pos_str = pos_str.replace(pos_str.split(
                        '-')[-1].split(':')[0], str(int(pos_str.split('-')[-1].split(':')[0])+1))
                elif pos_str.split(':')[-1] == '-':
                    modified_pos_str = pos_str.replace(re.findall(
                        r'\d+', pos_str)[0], str(int(re.findall(
                            r'\d+', pos_str)[0])-1))

            if loop_count <= maxnumber:
                loop_look(direction_flag, infasta, modified_pos_str, trans_flag,
                          loop_count, maxnumber, n, nuc_file_name, pro_file_name)
        # ################################################################################################################################
        # 第一层if else
        # 0为正确,1为起始x,2为内部错,3为末尾错
        elif seq_check_flag == 2:
            # pos_str
            # cds_seq, sorted_pos_list, flag_gene_type, len_trna_type
            # seq_check_flag, first_stop_codon_index_in_protein_seq, protein_seq
            # trans2protein_seq(cds_seq)函数对条件限定过了，因此以下代码条件比较宽松
            '''
            5070项目rpl16
            假设有俩终止子
            假设没有内含子
            假设第二段比第一段长
            '''
            if protein_seq.count('*') == 2:
                len_protein_seq_1, len_protein_seq_2 = get_current_first_end_pos(
                    sorted_pos_list, first_stop_codon_index_in_protein_seq)
                if (len_protein_seq_2 > len_protein_seq_1) and len(sorted_pos_list) == 1:
                    if pos_str.split(':')[-1] == '+':
                        modified_pos_str = pos_str.replace(pos_str.split(
                            '-')[0], str(int(pos_str.split('-')[0])+len_protein_seq_1))
                    elif pos_str.split(':')[-1] == '-':
                        modified_pos_str = pos_str.replace(re.findall(
                            r'\d+', pos_str)[-1], str(int(re.findall(
                                r'\d+', pos_str)[-1])-len_protein_seq_1))
                    if loop_count <= maxnumber:
                        loop_look(direction_flag,  infasta, modified_pos_str, trans_flag,
                                  loop_count, maxnumber, nuc_file_name, pro_file_name)
                '''
                5070项目ndhI
                假设有1终止子
                假设没有内含子
                第一段比第二段长
                '''
            elif protein_seq.count('*') == 1:
                len_protein_seq_1, len_protein_seq_2 = get_current_first_end_pos(
                    sorted_pos_list, first_stop_codon_index_in_protein_seq)
                if (len_protein_seq_2 < len_protein_seq_1) and len(sorted_pos_list) == 1 and maxnumber != 0:
                    if pos_str.split(':')[-1] == '+':  # 改末尾
                        modified_pos_str = pos_str.replace(pos_str.split(
                            '-')[-1], str(int(pos_str.split('-')[-1])-len_protein_seq_2))
                    elif pos_str.split(':')[-1] == '-':  # 改开头
                        modified_pos_str = pos_str.replace(re.findall(
                            r'\d+', pos_str)[0], str(int(re.findall(
                                r'\d+', pos_str)[0])+len_protein_seq_2))
                    if loop_count <= maxnumber:
                        loop_look(direction_flag,  infasta, modified_pos_str, trans_flag,
                                  loop_count, maxnumber, nuc_file_name, pro_file_name)
            else:
                print(sorted_pos_list)
                print(modified_pos_str)
                print('skip')
        # #########################################################################################################################################
        # 第一层if else
        # 0为正确,1为起始x,2为内部错,3为末尾错
        elif seq_check_flag == 0:  # 起始ok

            if len(pos_str.split(';')) != len(sorted_pos_list):  # ncbi上跨首尾基因会写成一段，实则两段
                if len(sorted_pos_list) == 2:
                    modified_pos_str = sorted_pos_list[0].replace(
                        ':1', ':+')+';'+sorted_pos_list[1].replace(':1', ':+')  # 考虑跨首尾
                else:
                    modified_pos_str = sorted_pos_list
            else:
                modified_pos_str = pos_str
            print('Correct Position: [{}]'.format(modified_pos_str))
            if pro_file_name != 'NULL':
                with open(os.path.join(current_abs_path, pro_file_name+'.protein_seq'), 'w') as f_handle:
                    f_handle.write(str(protein_seq)+'\n')
    print('[{}]'.format(modified_pos_str))
    return sorted_pos_list, first_stop_codon_index_in_protein_seq, modified_pos_str, seq_check_flag


###################################################################################################################


# 如果内部有终止子,则开始尝试返回新的基因位置，指开头到第一个终止子这一段


def get_current_first_end_pos(sorted_pos_list, first_stop_codon_index_in_protein_seq):
    # sorted_pos_list = []  # 排序后位置,起始子序列在列表第一位，终止子在列表最后一位
    # first_stop_codon_index_in_protein_seq = 200  # 包括第一个终止子在内的前面所有密码子个数
    # 包括第一个终止子在内的前面所有碱基数
    inter_pos = 3*(first_stop_codon_index_in_protein_seq+1)

    strand_list = []
    lenth_list = []
    for ele in sorted_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        lenth = end-start+1
        lenth_list.append(lenth)
        strand_list.append(strand)

    print('Lenth list:{0}\tStrand list:{1}\tTotal length:{2}bp'.format(
        lenth_list, strand_list, sum(i for i in lenth_list)))  # 序列可能是多段

    if inter_pos == sum(i for i in lenth_list):  # 等于全长说明序列正确，终止子在末尾
        print('Stop codon lie in [{}]'.format(
            sorted_pos_list[-1]))  # 说明终止子位于最后一段内
    else:
        remaining_bp = sum(i for i in lenth_list)-inter_pos  # 剩余的碱基数,不包括第一个终止子
        if remaining_bp <= lenth_list[-1]:  # 剩余的长度小于最后一段，说明终止子位于最后一段内
            print('Stop codon lie in [{}]'.format(sorted_pos_list[-1]))
            # 20220808 更正第一个终止密码子出现位置的计算公式
            if sorted_pos_list[-1].split(':')[-1] == 1:  # plus链基因
                current_first_end_pos = int(re.findall(
                    r'\d+', sorted_pos_list[-1])[-1])-remaining_bp-1  # 终止子中间那个碱基位置
                print('{}-{}:+'.format(current_first_end_pos -
                      1, current_first_end_pos+1))
                print('\n')
            elif sorted_pos_list[-1].split(':')[-1] == -1:  # minus链基因
                current_first_end_pos = int(re.findall(
                    r'\d+', sorted_pos_list[-1])[0])+remaining_bp+1
                print('{}-{}:-'.format(current_first_end_pos -
                      1, current_first_end_pos+1))
                print('\n')

        elif inter_pos <= lenth_list[0]:  # 包括第一个终止在内的长度小于第一段，说明终止子位于第一段
            print('Stop codon lie in [{}]'.format(sorted_pos_list[0]))
        else:  # 这里是简化了，不在第一段和最后一段，就认为在中间那段 ，可能存在bug
            print('Stop codon lie in [{}]'.format(sorted_pos_list[1]))
    return inter_pos, remaining_bp  # 前半截，后半截

#################################################################################################################


if __name__ == '__main__':

    loop_count = 0  # 控制递归次数,在loop_look函数外部定义全局变量   递归的计数
    loop_count_flag1 = []  # 20221020 定义一个查找正确起始子次数的列表,作为提示信息向外输出
    loop_count_flag3 = []
    loop_count_flag = []

    sorted_pos_list, first_stop_codon_index_in_protein_seq, \
        modified_pos_str, seq_check_flag = loop_look(args.direction_flag, args.infasta, args.pos_str,
                                                     args.trans_flag, loop_count, args.maxnumber,
                                                     args.codonnumber, args.nuc_file_name, args.pro_file_name)
    if type(first_stop_codon_index_in_protein_seq) == type(1):
        get_current_first_end_pos(
            format_pos(modified_pos_str), first_stop_codon_index_in_protein_seq)
