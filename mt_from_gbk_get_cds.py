#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_from_gbk_get_cds.py
#         Author:   yujie
#    Description:   mt_from_gbk_get_cds.py
#        Version:   1.0
#           Time:   2022/03/09 15:21:51
#  Last Modified:   2022/03/09 15:21:51
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
from Bio import SeqIO
import os
import re

parser = argparse.ArgumentParser(add_help=False, usage='\npython3   将fa序列反向互补')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[xxx.fasta]', help='输入fa文件', type=str, required=False)
optional.add_argument(
    '-l', '--lenth', metavar='[基因序列长度]', type=bool, help="有输入文件即可", default='1', required=False)


optional.add_argument('-s1', '--seq1',
                      metavar='[ATG→CAT]', help='DNA反向互补', type=str, required=False)
optional.add_argument('-s2', '--seq2',
                      metavar='[CAU→AUG]', help='RNA反向互补', type=str, required=False)
optional.add_argument('-s3', '--seq3',
                      metavar='[CAU→ATG]', help='RNA与DNA间反向互补', type=str, required=False)
optional.add_argument('-s4', '--seq4',
                      metavar='[U→T]', help='不反向不互补仅替换', type=str, required=False)
optional.add_argument('-s5', '--seq5',
                      metavar='[str.upper()]', help='字符串大写', type=str, required=False)
optional.add_argument('-s6', '--seq6',
                      metavar='[str.lower()]', help='字符串小写', type=str, required=False)

optional.add_argument('-o', '--output',
                      metavar='[ir_xxx.fasta]', help='输出反向后的fa文件', type=str, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note + format_seq + "\n"


def get_cds(gbk_file, f_cds):
    seq_record = SeqIO.read(gbk_file, "genbank")
    # print(seq_record.seq)
    complete_seq = str(seq_record.seq)
    # print(seq_record)
    # print(seq_record.annotations["accessions"])
    complete_note = ">" + seq_record.id + ":" + \
        seq_record.annotations["accessions"][0] + \
        " " + seq_record.description + "\n"
    complete_fasta = format_fasta(complete_note, complete_seq, 70)  # 70换行本例不采用
    # print(type(seq_record.features))

    cds_fasta = ""
    for ele in seq_record.features:
        if ele.type == "CDS":
            # print(ele.qualifiers)
            cds_seq = ""
            tmp_list = []
            # print(ele.location.parts)
            for ele1 in ele.location.parts:
                # print(ele1.start)
                # print(ele1.end)
                #print(int(re.findall(r'\d+', str(ele1.start))[0]))
                tmp_list.append(re.findall(r'\d+', str(ele1.start))[0])
                tmp_list.append(re.findall(r'\d+', str(ele1.end))[0])
                cds_seq += complete_seq[ele1.start:ele1.end]

            cds_note = ">" + seq_record.id + \
                " [" + str(int(tmp_list[0])+1)+".." + tmp_list[-1]+"]" + \
                " [gene=" + ele.qualifiers['gene'][0] + "]" + "\n"
            cds_fasta += format_fasta(cds_note, cds_seq, 70)
            print(cds_note)

            if (f_cds):  # ele有可能是trna,要确保先找到一个cds后才能退出,所以放上面if的下一级
                break
    return cds_fasta, complete_fasta


if __name__ == '__main__':
    # 文件输出路径
    out_cds_file_path = "F:/Hibiscus_sabdariffa/out/cds.fasta"
    out_complete_file = "F:/Hibiscus_sabdariffa/out/complete.fasta"
    # genbank 文件路径
    genbank_dir_path = "F:\\Hibiscus_sabdariffa\\111"
    out_cds_file_path_obj = open(out_cds_file_path, "w")
    out_complete_file_obj = open(out_complete_file, "w")
    for file in os.listdir(genbank_dir_path):
        # print(os.sep)
        cds_fasta, complete_fasta = get_cds(
            genbank_dir_path + os.sep + file, False)
        out_cds_file_path_obj.write(cds_fasta)
        out_complete_file_obj.write(complete_fasta)
