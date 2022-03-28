#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_blast_match.py
#         Author:   yujie
#    Description:   mt_blast_match.py
#        Version:   1.0
#           Time:   2022/03/11 09:49:54
#  Last Modified:   2022/03/11 09:49:54
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


createvar = locals()

parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
python3   mt_from_gbk_get_cds.py\n\
根据blast 结果\n\
输出infile每一行信息\n\
输出ref对应长度信息\n\
把提及到的序列放一起')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='sample-cds', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-r', '--ref',
                      metavar='[file]', help='ref-cds', type=str, default="F:\\ref_tre\\gene\\feature\\ref.gene.seq", required=False)
optional.add_argument('-o', '--outdir',
                      metavar='[dir]', help='输出的路径', type=str, default="F:\\ref_tre\\gene\\blast", required=False)
optional.add_argument('-f', '--flag',
                      metavar='[bool]', help='是否运行blast,默认否', type=bool, default=False, required=False)
optional.add_argument('-c', '--checkflag',
                      metavar='[bool]', help='是否设定阈值,默认否', type=bool, default=False, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

##################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))

if args.flag:
    cmd = "formatdb -i {} -p F -o F".format(args.infile)
    # print(cmd)
    os.system(cmd)

    cmd = "blastn -query {0} -db {1} -evalue 1e-5 -outfmt 5 -max_hsps 1 -out {2}/blastn.result.xml -word_size 7".format(
        args.ref, args.infile, args.outdir)
    # print(cmd)
    os.system(cmd)

    cmd = "perl /share/nas6/xul/program/mt2/phytree/gene_tree/src/blast_parser.pl -tophit 1 -topmatch 1 -m 7 {0}/blastn.result.xml > {0}/blastn.tophit.result.xls".format(
        args.outdir)
    # print(cmd)
    os.system(cmd)

"""
根据blast 结果
输出infile每一行信息
输出ref对应长度信息
把提及到的序列放一起
"""

#######################################################################################


def read_fasta_to_dic(infasta):
    with open(infasta, 'r') as f:
        seq_id = ''
        seq_index = []
        dict_seq = {}
        dict_len = {}
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip('\n')
                seq_index.append(line.replace(
                    "\n", "").replace(">", ""))
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{2} Item Quantity: {0} {1}'.format(
        len(dict_seq), len(dict_len), os.path.basename(infasta)))
    return dict_seq, dict_len


sample_cds = {}
sample_len = {}
sample_cds, sample_len = read_fasta_to_dic(args.infile)
ref_cds = {}
ref_len = {}
ref_cds, ref_len = read_fasta_to_dic(args.ref)  # >
for i in sample_cds.keys():
    createvar[i] = []
    # print(createvar[i])
print("Count Completed\n")

######################################################################################
blastn_tophit_result_path = os.path.join(
    args.outdir, 'blastn.tophit.result.xls')
cds_homo = {}

# 以二进制方式读取文件,统计换行符个数,即行数
count = 0
thefile = open(blastn_tophit_result_path, 'rb')
while True:
    buffer = thefile.read(1024 * 8192)
    if not buffer:
        break
    count += buffer.count('\n'.encode())
thefile.close()
print('line Total {}'.format(count))


with open(blastn_tophit_result_path, 'r') as f:
    for i in range(0, 1):
        f.readline()
    n = count-1  # n行有用信息
    count_n1 = 0
    count_n2 = 0
    line_number = 1
    for line in f:
        line_number += 1
        # print(line.strip())
        content = line.split('\t')
        (qid, q_length, q_aln_start, q_aln_end, s_length, s_aln_start,
         s_aln_end, sid) = content[0], int(content[1]), int(content[2]), int(content[3]), int(content[5]), int(content[6]), int(content[7]), content[15].strip()
        #cds_homo[">"+sid] = []
        #print('>{0} len: {1}'.format(sid.strip(), s_length))
        # subject_start_end_list = re.findall(r'\d+', content[15].split()[1])
        #print('{0} {1}'.format(qid.split()[0], q_length))
        # query_start_end_list = re.findall(r'\d+', content[0].split()[1])
        #print(q_aln_start > q_aln_end or s_aln_start > s_aln_end or (q_aln_end - q_aln_start + 1) / q_length < 0.6 or (s_aln_end - s_aln_start + 1) / s_length < 0.6)
        if q_aln_start > q_aln_end and args.checkflag:
            count_n2 += 1
            print("skip{0}行 {1}  {2}".format(
                line_number, qid, 'q_aln_start > q_aln_end'))
            next
        elif s_aln_start > s_aln_end and args.checkflag:
            count_n2 += 1
            print("skip{0}行 {1}  {2}".format(
                line_number, qid, 's_aln_start > s_aln_end'))
            next
        elif (q_aln_end - q_aln_start + 1) / q_length < 0.6 and args.checkflag:
            count_n2 += 1
            print("skip{0}行 {1}  {2}".format(line_number, qid,
                  '参考比自身的比率 < 0.6'))
            next
        elif (s_aln_end - s_aln_start + 1) / s_length < 0.6 and args.checkflag:
            count_n2 += 1
            print("skip{0}行 {1}  {2}".format(line_number, qid,
                  '(s_aln_end - s_aln_start + 1) / s_length < 0.6'))
            next
        else:
            # print(sid)
            # print(createvar['>'+sid])
            #print("当前为{}行 Efficient".format(line_number))
            count_n1 += 1
            q_pos = re.findall(r'\d+', qid.split()[1])  # 取位置出来
            q_start, q_end = q_pos[0], q_pos[1]
            #print(q_start, q_end)
            s_pos = re.findall(r'\d+', sid.split()[1])  # 取位置出来
            s_start, s_end = s_pos[0], s_pos[1]
            #print(s_start, s_end)
            dis = len(ref_cds[">"+qid])
            tmp = [qid.split()[0], dis, ">"+qid]  # ref信息
            # print(tmp)
            # print(createvar(sid))
            createvar['>'+sid].append(tmp)
            #print('{0}{1}'.format(sid, len(createvar['>'+sid])))
            cds_homo[">"+sid] = createvar['>'+sid]

print('Total {} gene'.format(len(cds_homo)))
if n == count_n1+count_n2:
    print('== {} {} {}'.format(n, count_n1, count_n2))
else:
    print("Please check!")
# print(len(cds_homo))
# print(cds_homo)
print('Parsing Completed\n')

#########################################################################################
output = []
homo_group = []
n = 0  # 最后写入了多少条序列
for k in cds_homo.keys():
    fasta = ''
    # print(k)
    NCBI_id, pos, gene = k.split()[0], re.findall(
        r'\d+', k.split()[1]), k.split()[2]
    #print(NCBI_id, pos, gene)
    start, end = pos[0], pos[1]
    gene = gene.split('=')[1].replace("]", '')
    # print(gene)
    fasta = (k+'\n'+sample_cds[k])
    # print(fasta)
    lenth = sample_len[k]  # G1长度
    my_filter = {}
    homo_group = cds_homo[k]
    # for i in homo_group:
    # tmp=
    homo_group.sort(key=lambda x: x[0])  # 以数组每个元素中的第二个元素排序
    # print(homo_group)  # .sort(key=takeSecond))
    print('\n')
    print('{0} len: {1} Total: {2}'.format(
        k, lenth, 1+len(homo_group)))  # 物种基因情况,total指包括物种1在内所有物种拥有该基因的物种个数
    for i in homo_group:
        print(i[0], i[1])
    n += 1+len(homo_group)
print(count_n1, len(sample_len), n)
print("Print Completed\n")
##########################################################################


def format_fasta(note, seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        # if (index + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return note + format_seq + "\n"


# sample_cds ref_cds  存有 >名字  对应 序列
# cds_homo 存有 MmG1 对应 其他参考的对应基因名字及长度
n = 0
for i in cds_homo.keys():
    n += 1
    gene = i.split()[2].split('=')[1].rstrip(']')
    filename = 'gene{0}.{1}.fasta'.format(n, gene)
    # print(filename)
    # print(cds_homo[i])
    with open(os.path.join(args.outdir, 'fasta', filename), 'w') as f:
        sample_fasta = i+'\n'+sample_cds[i]+'\n'
        f.write(sample_fasta)
        for j in cds_homo[i]:
            #fasta = format_fasta(j[2], ref_cds[j[2]], 70)
            ref_fasta = j[2]+'\n'+ref_cds[j[2]]+'\n'
            f.write(ref_fasta)

# print('\n')
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
###########
