#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   fungi_from_gbk_get_cds_V1.0.py
#         Author:   yujie
#    Description:   fungi_from_gbk_get_cds_V1.0.py
#        Version:   1.0
#           Time:   2023/02/01 14:02:07
#  Last Modified:   2023/02/01 14:02:07
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
# from icecream import ic
import argparse
# import linecache
import os
# import pretty_errors
import re
import sys
import time
# import copy
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
python3 cp_from_gbk_get_cds_V3.0.py -i [gbk dir] -o [out dir]\n\
Each species generates three files of cds,trna and complete sequence.')
optional = parser.add_argument_group('optional')
required = parser.add_argument_group('required')
required.add_argument(
    '-i', '--input',  help='gbk dir path', type=str)
required.add_argument(
    '-o', '--output',  help='output dir path', type=str)
optional.add_argument('-c', '--check', help='create genes via gap,not running by default,input "-c" when using',
                      action='store_true', required=False)
optional.add_argument('-d', '--duplicates', help='remove duplicates,default true',
                      action='store_false', required=False)
optional.add_argument('-ignore', '--ignore', help='ignore prompt,default true',
                      action='store_false', required=False)
optional.add_argument('-info', '--info', help='show update log and exit',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help',
                      help='show this help message and exit')
args = parser.parse_args()

if args.info:
    print('\n更新日志:')
    print('\t20220906 修改基因组的类型判断,若不符合则退出')
    print('\t20221012 修改基因名字映射函数,使其能够识别CO1样式')
    print('\t20221111 增加gbk里基因名是全称(如ATPase)时的处理')
    print('\t20221217 feat: ✨ 对gbk文件进行去重')
    print('\t20221219 🐞fix(get_gene_note): 修改没有gene标签的cds的默认ID')
    print('\t20221220 ✨feat(main): 写入去重后的登录号')
    print('\n')
    sys.exit(0)

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
# ###########################################################################################################


def get_complete_note(seq_record):  # 获取整个完整基因组ID
    try:
        seq_id = ''
        seq_record.description = seq_record.description.replace('/', '_')
        # 20220819 NC_044756.1.gbk voucher Liu HM/CP02 chloroplast  第一次处理,去掉特殊符号
        # if seq_record.description.find('chloroplast'):#有bug,用str格式化后就没问题了
        # 20220627 if str(seq_record.description).find('chloroplast') -1也成立,判断时一定要以False True为准
        # or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' or seq_record.description.split(',')[-2].split()[-1] == 'plastid':

        if str(seq_record.description).find('chloroplast') > 0 \
            or str(seq_record.description).find('plastid') > 0 \
                or seq_record.description.split(',')[-2].split()[-1] == 'chloroplast' \
                or seq_record.description.split(',')[-2].split()[-1] == 'plastid':  # gbk 属于 叶绿体
            if str(seq_record.description).find('chloroplast') > 0:  # 20220819 根据不同关键词分割
                seq_id = seq_record.description.split(
                    'chloroplast')[0].replace(' ', '_').rstrip('_')
            elif str(seq_record.description).find('plastid') > 0:
                seq_id = seq_record.description.split(
                    'plastid')[0].replace(' ', '_').rstrip('_')
            name = seq_record.name
            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"  # chloroplast--叶绿体

        elif str(seq_record.description).find('mitochondrion') > 0 \
                or str(seq_record.description).find('mitochondrial') > 0:
            # 20220811 NC_031548.gbk 描述部分 DEFINITION  Oxynoemacheilus angorae mitochondrial DNA, complete genome.
            seq_id = seq_record.description.split(
                'mitochondrion')[0].replace(' ', '_').rstrip('_')  # 形如 UNVERIFIED:_Rumina_decollata
            if seq_id.startswith('UNVERIFIED:_'):  # 第二次处理,去掉 UNVERIFIED:_
                # seq_id = seq_id.lstrip('UNVERIFIED:_')  # 有bug ???????
                seq_id = seq_id.split(':')[1].strip('_')
            # Cerion_watlingense_voucher_USNM:1514170_MN904501 第三次处理,去掉冒号后内容
            if len(seq_id.split(':')) > 1:
                seq_id = seq_id.split(':')[0]
            name = seq_record.name  # 要么是登录号  要么是样本
            if seq_id == name:
                seq_id = seq_id
            elif seq_id != name:
                seq_id = seq_id+'_'+name
            complete_note = ">" + seq_id + "\n"  # mitochondrion--线粒体

        else:
            print(seq_record.description)
            print('Genome Type WARNING! {}!\n{}'.format(
                seq_record.description.split(', ')[-2].split()[-1],
                'Program will automatically process!'))
            complete_note = ">" + (seq_record.description.split('mitochondrion')
                                   [0]).replace(' ', '_').rstrip('_') + "\n"

    except:  # 如果遇到任何出错
        print('try/except')
        complete_note = ''
        #gbk_type = input('genome type(1:chloroplast;2:mitochondrion): ')
        gbk_type = 2
        if gbk_type == 1:
            seq_id = seq_record.description.split(
                'chloroplast')[0].replace(' ', '_').rstrip('_')  # 物种或样品名
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
        elif gbk_type == 2:
            print(seq_record.description)
            seq_id = input('seq_id: ')
            complete_note = ">" + seq_id + "\n"

    return complete_note, seq_id

# ###################################################################################################################
# 仅在 get_gene_note(ele, complete_seq, seq_id, tmp_gene_name)中使用，
# 以应对"gbk中cds没有/gene标签，但是有/product标签"的情况


def gene_name_standardization_1(gene_name):  # 线粒体 格式化基因名字,可重复使用
    all_gene_dict = {'ATP synthase F0 subunit 6': 'ATP6', 'ATP synthase F0 subunit 8': 'ATP8', 'cytochrome b': 'CYTB',
                     'cytochrome c oxidase subunit I': 'COX1', 'cytochrome c oxidase subunit II': 'COX2',
                     'cytochrome c oxidase subunit III': 'COX3', 'NADH dehydrogenase subunit 1': 'ND1',
                     'NADH dehydrogenase subunit 2': 'ND2', 'NADH dehydrogenase subunit 3': 'ND3',
                     'NADH dehydrogenase subunit 4': 'ND4', 'NADH dehydrogenase subunit 4L': 'ND4L',
                     'NADH dehydrogenase subunit 5': 'ND5', 'NADH dehydrogenase subunit 6': 'ND6'}
    if gene_name in all_gene_dict.keys():
        gene_name = all_gene_dict[gene_name]
    return gene_name


def gene_name_standardization_2(gene_name):  # 格式化基因名字,可重复使用
    # 初始化
    all_gene_list_upper = ['ATP6', 'ATP8', 'CYTB', 'COX1', 'COX2',
                           'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_upper2 = ['ATP6', 'ATP8', 'CYTB', 'CO1', 'CO2',
                            'CO3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']
    all_gene_list_lower = ['atp6', 'atp8', 'cob', 'cox1', 'cox2',
                           'cox3', 'nad1', 'nad2', 'nad3', 'nad4', 'nad4l', 'nad5', 'nad6']
    name_flag = 0
    # 基因名字前处理1 # 20220825 gbk文件里 基因名是罗马数字，因此需要先处理
    gene_name = gene_name.replace('III', '3').replace(
        'II', '2').replace('I', '1')
    # 基因名字前处理2 # 20221111 gbk文件里 基因名是全称 如'ATPase '
    gene_name = gene_name.replace('ATPase ', 'ATP').replace('Cyt ', 'CYT')

    # ------------------名字映射
    # 名字样式大写的情况
    if gene_name.upper() in all_gene_list_upper:
        gene_name = gene_name.upper()
    else:
        i = 0
        while i < 13:
            # 20221012  MN053900.1.gbk 有一个CO1
            if all_gene_list_upper2[i] == gene_name.upper():
                gene_name = all_gene_list_upper[i]
                print(gene_name)
                break
            else:
                i += 1
        if i >= 13:
            #print(gene_name, 'Try further changes...')
            # 名字样式小写的情况
            i = 0
            while i < 13:
                # 20220624  JN619347.gbk 有一个nad4L 大小写混合 就离谱
                if all_gene_list_lower[i] == gene_name.lower():
                    gene_name = all_gene_list_upper[i]
                    break
                else:
                    i += 1
            if i >= 13:
                if args.ignore == False:
                    print('all_three_methods_failed!')
                    print(gene_name, 'WARNING!Please check!')
                name_flag = 1
    return gene_name, name_flag

# ###################################################################################################################


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
    return tmp_list, cds_seq


# ###################################################################################################################
def get_gene_note(file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name):  # 获取gene的id及序列
    '''
    传入上一个基因信息tmp_gene_name,返回当前基因信息tmp_gene_name
    OrderedDict(
        [
        ('codon_start', ['1']),
        ('gene', ['rps19']),
        ('product', ['ribosomal protein S19']),
        ('protein_id', ['UKP82335.1']),
        ('transl_table', ['11']),
        ('translation', [
         'MTRSLKKNPFVANHLLRKINKLNTKAEKEIIITWSRASTIIPTMIGHTIAIHNGKEHLPIYITDRMVGHKLGEFSPTLNFRGHAKNDNRSRR'])
        ]
        )
    '''
    # 20220825 NC_034226.gbk cds没有/gene标签，但是有/product标签
    if 'gene' not in ele.qualifiers.keys():  # 返回上一个基因,好从其他参考找这个没名字的
        if 'product' in ele.qualifiers.keys():
            tmp_gene_name = ele.qualifiers['product'][0]
            tmp_gene_name = gene_name_standardization_1(tmp_gene_name)
        else:
            try:
                tmp_gene_name = tmp_gene_name+'_next'
            except:
                tmp_gene_name = input(
                    "\n{0}: {1} Previous: {2}. Current: {3}.\nPlease input current gene name:".
                    format(file_no, file_name, tmp_gene_name, ele.location.parts))

    elif 'gene' in ele.qualifiers.keys():
        tmp_gene_name = ele.qualifiers['gene'][0]

    tmp_gene_name, name_flag = gene_name_standardization_2(
        tmp_gene_name)  # 20221012 也需要第二种标准化

    tmp_list, gene_seq = merge_sequence(ele, complete_seq)
    pos_info_in_gene_note = " ["
    for i in range(2*len(ele.location.parts)):
        if i % 2 == 0:
            pos_info_in_gene_note += tmp_list[i]+'..'
        elif i % 2 == 1:
            pos_info_in_gene_note += tmp_list[i]+';'
    pos_info_in_gene_note = pos_info_in_gene_note.rstrip(';')+"]"
    gene_note = ">" + seq_id + pos_info_in_gene_note + " [gene=" + \
        tmp_gene_name + "]" + "\n"  # '>'后的格式和已有脚本兼容
    return gene_note, gene_seq, tmp_gene_name

# ##############################################################################################


def get_gene(gbk_file_path, flag, dict_gene_len, file_no):  # 解析gbk文件获取cds
    """完整基因组"""
    file_name = os.path.basename(gbk_file_path)
    seq_record = SeqIO.read(gbk_file_path, "genbank")
    complete_seq = str(seq_record.seq)
    complete_note, seq_id = get_complete_note(seq_record)
    complete_fasta = format_fasta(complete_note, complete_seq, 70)  # 70换行本例不采用
    """gene序列"""
    cds_count = 0  # 对cds数量计数
    trna_count = 0
    cds_fasta = ""
    trna_fasta = ""
    list_cds_name = []  # 统计cds种类，记录重复
    list_cds_name_len = []  # 统计cds种类，记录重复,记录长度
    list_trna_name = []  # 统计trna种类,记录重复
    tmp_gene_name = ''  # 上一个基因名字,为子函数get_gene_note()准备的
    for ele in seq_record.features:
        '''
        CDS
        tRNA
        rRNA
        repeat_region
        '''
        if ele.type == "CDS":
            cds_count += 1
            cds_note, cds_seq, tmp_gene_name = get_gene_note(
                file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name)
            # list_cds_name.append(tmp_gene_name)  # 本次的基因名字 复用,线粒体的话,在下一部分存入列表
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            gene_name = tmp_gene_name
            gene_name, name_flag = gene_name_standardization_2(gene_name)
            list_cds_name.append(gene_name)  # 存入列表
            list_cds_name_len.append(
                gene_name+':'+str(len(cds_seq)))  # 存入列表
            if name_flag == 1:
                print(gbk_file_path, 'WARNING!Please check!')
            dict_gene_len[gene_name].append(
                3*(len(ele.qualifiers['translation'][0])+1))  # cds序列长度
            if (flag):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
        elif ele.type == 'tRNA':
            trna_count += 1
            trna_note, trna_seq, tmp_gene_name = get_gene_note(
                file_no, file_name, ele, complete_seq, seq_id, tmp_gene_name)
            trna_fasta += format_fasta(trna_note, trna_seq, 70)
            gene_name = tmp_gene_name
            gene_name = gene_name  # gene_name_standardization(gene_name)
            list_trna_name.append(gene_name)  # 存入列表

    s = '{2}: {0} has {1} CDS'.format(file_name, cds_count, file_no)
    if cds_count == 0:
        # --------There may be no comments--------'.format(
        s = '{2}: {0} has {1} CDS'.format(file_name, cds_count, file_no)
        print(s.ljust(50), '----------There may be no comments----------')
    elif cds_count != 0:
        # 20220811 输出左对齐  str.ljust(50)  达到50个字符宽度
        if args.ignore == False:
            print(s.ljust(50), '+', list_cds_name_len)
        else:
            print(s.ljust(50), '+', list_cds_name_len)
    return file_name, seq_id, complete_fasta, cds_fasta, cds_count, list_cds_name,  trna_fasta, trna_count, list_trna_name, dict_gene_len, s, list_cds_name_len

# ###############################################################################构造子函数


def create_gene_by_gap(dict_missing_gene, dict_gene_len, out_cds_file_list):  # 用gap构造没有的基因
    # dict_missing_gene   '>Achatina_fulica_NC_024601': []
    # out_cds_file_list   'cds_NC_036381.fasta'
    out_cds_file_list_1 = []
    for i in dict_missing_gene.keys():
        cds_fasta = ''
        for j in dict_missing_gene[i]:
            ave = round(sum(dict_gene_len[j]) /
                        len(dict_gene_len[j]))  # 该基因平均长度
            cds_note = (i+' [0..0]'+' [gene={}]').format(j)
            cds_seq = ave*'-'
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
        for m in out_cds_file_list:
            m = m.lstrip('cds_').rstrip('.fasta')
            if i.find(m) >= 0 or i.find(m.rstrip('.1')) >= 0:  # i里头也可能不带版本号,去掉m的再查找
                cds_file_path = 'cds_'+m+'.fasta'
                out_cds_file_list_1.append(cds_file_path)

        with open(args.output+os.sep+cds_file_path, 'ab+') as f_cds:  # 读写打开一个二进制文件，允许读，或在文件末追加数据
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
        os.makedirs(args.output)
    if (not os.path.exists(os.path.join(args.output, 'trna'))):
        os.makedirs(os.path.join(args.output, 'trna'))
    if (not os.path.exists(os.path.join(args.output, 'cds'))):
        os.makedirs(os.path.join(args.output, 'cds'))
    if (not os.path.exists(os.path.join(args.output, 'complete'))):
        os.makedirs(os.path.join(args.output, 'complete'))

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
    list_seq_id_de_duplication = []  # 去重的物种
    list_unique_accession = []  # 去重的登录号
    for i in all_gene_list_upper:
        dict_gene_len[i] = []

    """初始化"""
    dict_file_cds_count = {}  # 每个文件中cds计数
    file_list = [x for x in os.listdir(
        args.input) if (os.path.isfile(os.path.join(args.input, x)) and x.endswith('gbk'))]
    file_list.sort()  # key=lambda x: int(x.split('.')[0])) #根据文件名中的数字

    """主程序"""
    file_no = 0  # 为了屏幕输出更直观,在文件前面的编号
    out_cds_file_list = []  # 所有输出的cds前缀文件路径,用于用GAP构造序列
    for file in file_list:
        file_no += 1
        gbk_file_path = os.path.join(args.input, file)
        file_name, seq_id, complete_fasta, cds_fasta, cds_count, list_cds_name,  \
            trna_fasta, trna_count, list_trna_name, dict_gene_len, s, list_cds_name_len = get_gene(
                gbk_file_path, False, dict_gene_len, file_no)
        dict_file_cds_count[file_name] = cds_count  # 每个文件中cds计数
        '''
        20221217对gbk文件去重
        '''
        seq_id_content = seq_id.split('_')
        if seq_id.find('NC_') > 0:
            species = '_'.join(seq_id_content[:-2])
            accession = '_'.join(seq_id_content[-2:])
        else:
            species = '_'.join(seq_id_content[:-1])
            accession = '_'.join(seq_id_content[-1:])
        if args.duplicates == True:  # 执行去重
            if species not in list_seq_id_de_duplication:
                list_seq_id_de_duplication.append(species)
                list_unique_accession.append(accession)

                """写入文件"""
                out_cds_file = ('cds_' + file_name.rstrip('.gbk')+'.fasta')
                out_cds_file_list.append(out_cds_file)  # gap构造基因
                out_complete_file = (seq_id+'.fasta')
                out_trna_file = ('trna_'+file_name.rstrip('.gbk')+'.fasta')
                out_log_file = 'log'
                with open((os.path.join(args.output, 'complete', out_complete_file)), 'wb') as f_complete, \
                    open((os.path.join(args.output, 'cds', out_cds_file)), 'wb') as f_cds,\
                        open(os.path.join(args.output, 'trna', out_trna_file), 'wb') as f_trna,\
                        open(os.path.join(args.output,  out_log_file), 'a+') as f_log:
                    f_complete.write(complete_fasta.encode())
                    f_cds.write(cds_fasta.encode())
                    f_trna.write(trna_fasta.encode())

                    """以下为统计部分"""
                    f_log.write(s+'\n')
                    f_log.write('>'+file.rstrip('.gbk')+'\t')
                    list_missing_gene = []  # 每个文件中缺失的基因统计,单列表
                    for i in range(len(all_gene_list_upper)):
                        if all_gene_list_upper[i] in list_cds_name \
                                or all_gene_list_upper[i].lower() in list_cds_name \
                                or all_gene_list_lower[i] in list_cds_name \
                                or all_gene_list_lower[i].upper() in list_cds_name:
                            f_log.write(all_gene_list_upper[i]+'\t')
                        else:
                            f_log.write('NULL'+'\t')
                            list_missing_gene.append(
                                all_gene_list_upper[i])  # 缺失的基因
                    if 1 <= len(list_missing_gene) < 13:
                        print(''.ljust(55), '-', list_missing_gene)

                    f_log.write('\n')
                    [f_log.write(tmp+'\t') for tmp in list_missing_gene]
                    f_log.write('\n')
                dict_missing_gene['>'+seq_id] = list_missing_gene  # 放入字典当中

        else:

            """写入文件"""
            out_cds_file = ('cds_' + file_name.rstrip('.gbk')+'.fasta')
            out_cds_file_list.append(out_cds_file)  # gap构造基因
            out_complete_file = (seq_id+'.fasta')
            out_trna_file = ('trna_'+file_name.rstrip('.gbk')+'.fasta')
            out_log_file = 'log'
            with open((os.path.join(args.output, 'complete', out_complete_file)), 'wb') as f_complete, \
                open((os.path.join(args.output, 'cds', out_cds_file)), 'wb') as f_cds,\
                    open(os.path.join(args.output, 'trna', out_trna_file), 'wb') as f_trna,\
                    open(os.path.join(args.output,  out_log_file), 'a+') as f_log:
                f_complete.write(complete_fasta.encode())
                f_cds.write(cds_fasta.encode())
                f_trna.write(trna_fasta.encode())

                """以下为统计部分"""
                f_log.write(s+'\n')
                f_log.write('>'+file.rstrip('.gbk')+'\t')
                list_missing_gene = []  # 每个文件中缺失的基因统计,单列表
                for i in range(len(all_gene_list_upper)):
                    if all_gene_list_upper[i] in list_cds_name \
                            or all_gene_list_upper[i].lower() in list_cds_name \
                            or all_gene_list_lower[i] in list_cds_name \
                            or all_gene_list_lower[i].upper() in list_cds_name:
                        f_log.write(all_gene_list_upper[i]+'\t')
                    else:
                        f_log.write('NULL'+'\t')
                        list_missing_gene.append(
                            all_gene_list_upper[i])  # 缺失的基因
                if 1 <= len(list_missing_gene) < 13:
                    print(''.ljust(55), '-', list_missing_gene)

                f_log.write('\n')
                [f_log.write(tmp+'\t') for tmp in list_missing_gene]
                f_log.write('\n')
            dict_missing_gene['>'+seq_id] = list_missing_gene  # 放入字典当中

    # 后面考虑在这给字典排个序
    with open((args.output+os.sep+'log'), 'a+') as f_log:
        f_log.write(str(dict_missing_gene))
    """除了物种1外,其他物种所有cds总数"""
    total_ref_gene = 0
    for i in dict_file_cds_count.keys():  # 键为seq_id,值为个数
        total_ref_gene += dict_file_cds_count[i]
    print('\n')
    # print(dict_gene_len)  # 键为每个基因,值为列表,列表为每个基因在不同物种中的长度
    # print(dict_missing_gene)  # 键为>seq_id,值为列表,列表为每个物种确实的基因
    # print(2*(total_ref_gene-13))
    # ic(len(dict_missing_gene))
    # ic(dict_missing_gene)
    # ic(len(out_cds_file_list))
    """gap构造基因"""
    if args.check:
        create_gene_by_gap(dict_missing_gene, dict_gene_len, out_cds_file_list)

    if len(list_unique_accession) != 0:
        print('\n{} left after removing duplicates\n'.format(
            len(list_unique_accession)))
    else:
        print('\nall species are reserved\n')

    with open((os.path.join(args.output, 'list_tre')), 'w') as list_handle:
        for i in list_unique_accession:
            list_handle.write(i+'\n')
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
