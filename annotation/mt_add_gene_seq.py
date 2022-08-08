#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_add_gene_seq.py
#         Author:   yujie
#    Description:   mt_add_gene_seq.py
#        Version:   2.0
#           Time:   2022/05/23 16:35:19
#  Last Modified:   2022/05/23 16:35:19
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import pretty_errors
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import re
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   mt_add_gene_seq.py\n\
\n\
1.常规使用\n\
1.1查看密码子 -n -i -p \n\
\n\
2.递归查找与存储\n\
2.1起始子查找,-m 最大查找次数\n\
2.2存储序列,-sn 基因名\n\
2.3存储蛋白,-sp 基因名\n\
\n\
Path: E:\OneDrive\jshy信息部\Script\mitochondrion\annotation\mt_add_gene_seq.py\n\
Version: V2.0'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infasta', metavar='[infasta]', help='输入fasta文件', type=str, default='E:\\Examples\\mt_from_gbk_get_cds\\Meghimatium_bilineatum_1_FULLMT.fsa', required=False)
optional.add_argument(
    '-p', '--posstr', metavar='[pos_str]', help="输入位置,形如'124353-124892:-;126001-126552:-'", type=str, default='14323-1527:+', required=False)
optional.add_argument(
    '-n', '--codonnumber', metavar='[codon_number]', help='密码子表,默认5', type=int, default=5, required=False)
optional.add_argument(
    '-m', '--maxnumber', metavar='[max_number]', help='最大递归查找次数,默认0,假查找', type=int, default=0, required=False)
optional.add_argument('-trans', '--trans_flag',
                      help='翻译?默认是,不运行则-c1', action='store_false', required=False)
optional.add_argument('-sn', '--nuc_file_name',
                      metavar='[store 2 dna]', help='默认否,值为NULL,存储则输入gene名', type=str,  default='NULL', required=False)
optional.add_argument('-sp', '--pro_file_name',
                      metavar='[store 2 protein]', help='默认否,值为NULL,存储则输入蛋白名', type=str,  default='NULL', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

##########################################################################################
# 子函数 功能简单


def read_file(infasta):  # 读取文件
    with open(infasta, 'r') as f:
        seq = ''
        for line in f:
            if not line.startswith('>'):
                seq += line.strip('\n')
    return seq


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


def format_pos(pos_str):  # 读取输入的位置为位置列表
    pos_list = []
    content = pos_str.split(';')
    for ele in content:
        if ele.split(':')[-1] == '-':
            tmp = ele.split(':')[0]+':'+'-1'
            pos_list.append(tmp)
        elif ele.split(':')[-1] == '+':
            tmp = ele.split(':')[0]+':'+'1'
            pos_list.append(tmp)
    return pos_list

#######################################################################################################################


def merge_sequence(pos_list, seq):  # 合并获取到的序列,顺便排一下位置顺序
    # pos_list 某基因格式化的位置
    # seq 全长序列

    # 20220629新增
    # 判断是否是trna,返回一个flag
    flag_gene_type = 'NULL'
    len_trna_type = 0
    if len(pos_list) == 1:
        start = pos_list[0].split(':')[0].split('-')[0]
        end = pos_list[0].split(':')[0].split('-')[-1]
        len_trna_type = abs(int(end)-int(start))+1
        if 55 <= len_trna_type <= 100:
            flag_gene_type = 'trna'

    # -----------------20220523 解决跨首尾基因
    seq_len = len(seq)
    if int(pos_list[0].split(':')[-1]) == 1 and int(pos_list[0].split(':')[0].split('-')[0]) > int(pos_list[0].split(':')[0].split('-')[-1]):  # 14323-1527:1
        pos1 = '{0}-{1}:1'.format(pos_list[0].split(':')
                                  [0].split('-')[0], seq_len)
        pos2 = '1-{}:1'.format(pos_list[0].split(':')[0].split('-')[-1])
        pos_list = [pos1, pos2]
    # ------------------
    cds_seq = ""
    tmp_pos_list = pos_list
    if int(pos_list[0].split(':')[-1]) == -1:  # 一般来说,有内含子的基因,几段方向相同,因此只判断第一段,完了重新排序即可
        tmp_pos_list = pos_list[::-1]  # 此时 起始子序列在列表第一位，终止子在最后一位

    for ele in tmp_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        start_index = start-1
        end_index = end
        if strand == (-1):
            # seq[start_index:end_index] 角标从start_index到end_index    取的是索引start-1一直到end  取的是start一直到end的碱基
            cds_seq += ir(seq[start_index:end_index])
        elif strand == (1):
            cds_seq += seq[start_index:end_index]
    return cds_seq, tmp_pos_list, flag_gene_type, len_trna_type

#######################################################################################################################


def trans2acid(cds_seq, n):  # 翻译成氨基酸,返回是否正确以及第一个终止子在基因序列上的相对位置
    # start_codon_list = ['TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'] #11

    # 20220601
    if n == 5:
        start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 5
    elif n == 2:
        start_codon_list = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
        end_codon_list = ['TAA', 'TAG', 'AGA',
                          'AGG', 'TA', 'T', 'AG']  # 2

    tmp_flag = 0  # tmp_flag 起始是否正确的标志,默认False   20220610改为数字,0为正确,1为起始x,2为内部错,3为末尾错
    inter_number = 0
    if len(cds_seq) % 3 == 1:
        print('len(sequence) not a multiple of three! {}=3n+1'.format(len(cds_seq)))
    elif len(cds_seq) % 3 == 2:
        print('len(sequence) not a multiple of three! {}=3n+2'.format(len(cds_seq)))

    coding_dna = Seq(cds_seq)
    acid = coding_dna.translate(table=n)
    print('------------------------------------------------------------')
    print(acid)

    if not cds_seq[0:3] in start_codon_list:  # 起始不正确
        print('#####start is wrong!')
        tmp_flag = 1
    else:  # 起始正确,又分3种情况
        if acid.count('*') > 1:  # 终止多于1,意味着提前终止
            print('#####interior is wrong!')
            tmp_flag = 2
            inter_number = acid.find('*')
            print(inter_number)
            print('\n')
        elif acid.count('*') < 1:  # 终止小于1 1.真的未终止2.线粒体终止了,共5种细分情况
            if len(cds_seq) % 3 == 1 and cds_seq[-1] in end_codon_list:
                print(
                    '------------------------------------------------------------ok')
            elif len(cds_seq) % 3 == 2 and cds_seq[-2:] in end_codon_list:
                print(
                    '------------------------------------------------------------ok')
            else:
                print('#####end is wrong!')
                tmp_flag = 3
        else:  # 1个终止子,提前终止或者没问题
            if not acid.endswith('*'):
                print('#####interior is wrong!')
                tmp_flag = 2
                inter_number = acid.find('*')
                print(inter_number)
                print('\n')
            else:
                tmp_flag = 0
                print('------------------------------------------------------------ok')
    return tmp_flag, inter_number, acid


###################################################################################################################


# 如果内部有终止子,则开始尝试返回新的基因位置，指开头到第一个终止子这一段
def get_current_first_end_pos(tmp_pos_list, inter_number):
    # tmp_pos_list = []  # 排序后位置,起始子序列在列表第一位，终止子在列表最后一位
    # inter_number = 200  # 包括第一个终止子在内的前面所有密码子个数
    inter_pos = 3*inter_number  # 包括第一个终止子在内的前面所有碱基数

    strand_list = []
    lenth_list = []
    for ele in tmp_pos_list:  # ele 1-10:-1
        strand = int(ele.split(':')[-1])
        start = int(ele.split(':')[0].split('-')[0])
        end = int(ele.split(':')[0].split('-')[-1])
        lenth = end-start+1
        lenth_list.append(lenth)
        strand_list.append(strand)

    print(lenth_list)  # 序列可能是多段
    print(strand_list)  # 序列可能是多段
    lenth_sum = 0  # 序列总长
    for i in lenth_list:
        lenth_sum += i
    remaining_bp = lenth_sum-inter_pos  # 剩余的碱基数,不包括第一个终止子

    if remaining_bp <= lenth_list[-1]:  # 剩余的长度小于最后一段，说明终止子位于最后一段内
        print('lie in [{}]'.format(tmp_pos_list[-1]))
        # 20220808 更正第一个终止密码子出现位置的计算公式
        if tmp_pos_list[-1].split(':')[-1] == 1:  # plus链基因
            current_first_end_pos = int(re.findall(
                r'\d+', tmp_pos_list[-1])[-1])-remaining_bp-1  # 终止子中间那个碱基位置
            print('{}-{}:+'.format(current_first_end_pos-1, current_first_end_pos+1))
            print('\n')
        elif tmp_pos_list[-1].split(':')[-1] == -1:  # minus链基因
            current_first_end_pos = int(re.findall(
                r'\d+', tmp_pos_list[-1])[0])+remaining_bp+1
            print('{}-{}:-'.format(current_first_end_pos-1, current_first_end_pos+1))
            print('\n')

    elif inter_pos <= lenth_list[0]:  # 包括第一个终止在内的长度小于第一段，说明终止子位于第一段
        print('lie in [{}]'.format(tmp_pos_list[0]))
    else:
        print('lie in [{}]'.format(tmp_pos_list[1]))
    return 0

#################################################################################################################

# 20220722 新增子函数


def storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq):  # 存储获取到的dna序列
    if flag_gene_type == 'trna':  # 20220629   trna 存起来
        print('\nType: tRNA  Len: '+str(len_trna_type)+'\n')
    # 20220808 进行精简
    current_abs_path = os.getcwd()
    if nuc_file_name != 'NULL':
        with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
            f_handle.write(cds_seq+'\n')

    # if flag_gene_type == 'NULL':  # 20220722   把 cds 存起来
        #current_abs_path = os.getcwd()
        # if nuc_file_name != 'NULL':
            # with open(os.path.join(current_abs_path, nuc_file_name), 'w') as f_handle:
            # f_handle.write(cds_seq+'\n')
    return 0

#################################################################################################################


# 循环查找   *.fas/"1-10:-;20-30:-"/翻译/递归计数/最大递归次数
def loop_look(infasta, posstr, trans_flag, loop_count, maxnumber, n, nuc_file_name, pro_file_name):
    inter_number = False  # 20220629 add  初始值为false

    seq = read_file(infasta)
    pos_list = format_pos(posstr)
    cds_seq, tmp_pos_list, flag_gene_type, len_trna_type = merge_sequence(
        pos_list, seq)  # tmp_pos_list  把位置当列表再传出来,这个位置信息向下传递
    print('\n'+cds_seq)
    print(tmp_pos_list)
    print(posstr)
    storage_dna(flag_gene_type, len_trna_type, nuc_file_name, cds_seq)

    if trans_flag and (flag_gene_type != 'trna'):  # 翻译
        tmp_flag, inter_number, acid = trans2acid(cds_seq, n)
        current_abs_path = os.getcwd()
        if pro_file_name != 'NULL':
            with open(os.path.join(current_abs_path, pro_file_name+'.acid'), 'w') as f_handle:
                f_handle.write(str(acid)+'\n')
        # #########################################################################################################################################
        # 第一层if else
        if tmp_flag == 0:  # 起始ok
            if len(posstr.split(';')) != len(tmp_pos_list):  # ??????
                if len(tmp_pos_list) == 2:
                    new_pos_str = tmp_pos_list[0].replace(
                        ':1', ':+')+';'+tmp_pos_list[1].replace(':1', ':+')  # 考虑跨首尾
                else:
                    new_pos_str = tmp_pos_list
            else:
                new_pos_str = posstr
            print('正确位置: {}'.format(new_pos_str))
        # #####################################################################################################################################
        # 第一层if else
            """考虑细分情况 20220610考虑起始子错误的查找  其他错误类型暂时不考虑,用原来的程序写死"""
        elif tmp_flag == 1:  # 起始错,这个优先要满足的条件就不对
            # posstr 能传到这里  形如 1-7:+;14020-14078:+
            # tmp_pos_list 也能传到这里  形如['1-7:1', '14020-14078:1']
            loop_count += 1
            print('第{}次查找中'.format(loop_count))
            cds_seq = cds_seq[3:]  # 已经判断起始错误了,因此直接把序列剪掉前面3个碱基

            if n == 5:
                start_codon_list = ['TTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG']
                end_codon_list = ['TAA', 'TAG', 'TA', 'T']  # 5
            elif n == 2:
                start_codon_list = ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']
                end_codon_list = ['TAA', 'TAG', 'AGA',
                                  'AGG', 'TA', 'T', 'AG']  # 2,转录时要加A
            # ##############################################################
            # 定义为第二层if else
            if cds_seq[0:3] not in start_codon_list and maxnumber != 0:  # 20220805  如果为假查找，就不进行下一步了
                start_flag = False
                print(tmp_pos_list)
                print(posstr)
                # 20220808 以下自动返回位置，也就是开头往后挪6bp
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[0], str(int(posstr.split('-')[0])+6))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[-1], str(int(re.findall(
                            r'\d+', posstr)[-1])-6))
                '''
                new_pos_str = input(
                    '与上次命令行输入-6bp new pos(末尾输.可重输): ')  # 先改手动输入,以后改自动
                while new_pos_str.find('.') > 0:  # 20220808
                    new_pos_str = input(
                        '与上次命令行输入-6bp new pos(末尾输.可重输): ')
                    if new_pos_str.find('.') < 0:
                        continue
                '''
                if loop_count <= maxnumber:
                    loop_look(infasta, new_pos_str, trans_flag,
                              loop_count, maxnumber, n, nuc_file_name, pro_file_name)
            elif cds_seq[0:3] in start_codon_list:
                start_flag = True
                print('\n'+cds_seq)
                # 20220808 以下自动返回正确位置，也就是开头往后挪3bp
                if posstr.split(':')[-1] == '+':
                    new_pos_str = posstr.replace(posstr.split(
                        '-')[0], str(int(posstr.split('-')[0])+3))
                elif posstr.split(':')[-1] == '-':
                    new_pos_str = posstr.replace(re.findall(
                        r'\d+', posstr)[-1], str(int(re.findall(
                            r'\d+', posstr)[-1])-3))
                tmp_flag, inter_number, acid = trans2acid(cds_seq, n)
                print('正确位置: {}'.format(new_pos_str))
        # ################################################################################################################################
        # 第一层if else
        else:
            maxnumber == 0  # 赋值为0  相当于一个假查找
            # 一般是在外面赋值,这里因为修改tmp_falg==1的查找,不能在参数设置时默认为0,否则会影响tmp_falg==1的情况
            loop_count += 1
            print('第{}次查找中'.format(loop_count))
            new_pos_str = '124353-124892:-;126001-126552:-'

            if loop_count <= maxnumber:
                loop_look(infasta, new_pos_str, trans_flag,
                          loop_count, maxnumber, n, nuc_file_name, pro_file_name)
            else:
                print('{}次查找未有结果,取消第{}次查找'.format(loop_count-1, loop_count))
    return tmp_pos_list, inter_number


if __name__ == '__main__':
    loop_count = 0  # 控制递归次数,在loop_look函数外部定义全局变量   递归的计数
    tmp_pos_list, inter_number = loop_look(
        args.infasta, args.posstr, args.trans_flag, loop_count, args.maxnumber, args.codonnumber, args.nuc_file_name, args.pro_file_name)
    if type(inter_number) == type(1):
        get_current_first_end_pos(tmp_pos_list, inter_number)

"""
def trans2acid(codon):  # 翻译成氨基酸
    """"""
    The Bacterial and Plant Plastid Code (11):
    Stnd    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    This    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts      = ---M---------------M------------MMMM---------------M------------
    Base1       = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2       = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3       = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """"""
    genetic_code_number = 11
    acid = ''
    code_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', }
    acid = code_table[codon]
    return acid
"""
