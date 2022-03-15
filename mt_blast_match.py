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

import argparse
from Bio import SeqIO
import os
import re
createvar = locals()

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   mt_from_gbk_get_cds.py')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='sample-cds', type=str, default="F:\\ref_tre\\gene\\feature\\Mm_G1.gene.seq", required=False)
optional.add_argument('-r', '--ref',
                      metavar='[file]', help='ref-cds', type=str, default="F:\\ref_tre\\gene\\feature\\ref.gene.seq", required=False)
optional.add_argument('-o', '--outdir',
                      metavar='[dir]', help='输出的路径', type=str, default="F:\\ref_tre\\gene\\blast", required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

cmd = "formatdb -i {} -p F -o F".format(args.infile)
# print(cmd)
# os.system(cmd)

cmd = "blastn -query {0} -db {1} -evalue 1e-5 -outfmt 5 -max_hsps 1 -out {2}/blastn.result.xml -word_size 7".format(
    args.ref, args.infile, args.outdir)
# print(cmd)
# os.system(cmd)

cmd = "perl /share/nas6/xul/program/mt2/phytree/gene_tree/src/blast_parser.pl -tophit 1 -topmatch 1 -m 7 {0}/blastn.result.xml > {0}/blastn.tophit.result.xls".format(
    args.outdir)
# print(cmd)
# os.system(cmd)


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
        print(len(dict_seq))
        print(len(dict_len))
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
with open(blastn_tophit_result_path, 'r') as f:
    for i in range(0, 1):
        f.readline()
    n = 260-1
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
        if(q_aln_start > q_aln_end or s_aln_start > s_aln_end or (q_aln_end - q_aln_start + 1) / q_length < 0.6 or (s_aln_end - s_aln_start + 1) / s_length < 0.6):
            count_n2 += 1
            # print("skip{}行".format(line_number))
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
    print(len(cds_homo))

    if n == count_n1+count_n2:
        print('==', n, count_n1, count_n2)
# print(len(cds_homo))
# print(cds_homo)
print('Parsing Completed\n')

#########################################################################################
output = []
homo_group = []
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
    print('{0} len: {1} Total: {2}'.format(
        k, lenth, 1+len(homo_group)))  # 物种基因情况
    for i in homo_group:
        print(i[0], i[1])
print("Print Completed")
