
import os
from Bio import SeqIO
import re


def format_fasta(ana, seq, num):
    """
    格式化文本为 fasta格式
    :param ana: 注释信息
    :param seq: 序列
    :param num: 序列换行时的字符个数
    :return: fasta格式文本
    """
    format_seq = ""
    for i, char in enumerate(seq):
        format_seq += char
        # if (i + 1) % num == 0:#可以用来换行
        #format_seq += "\n"
    return ana + format_seq + "\n"


def get_cds(gb_file, f_cds):
    """
    从 genbank 文件中提取 cds 序列及其完整序列
    :param gb_file: genbank文件路径
    :param f_cds: 是否只获取一个 CDS 序列
    :return: fasta 格式的 CDS 序列， fasta 格式的完整序列 
    """
    # 提取完整序列并格式为 fasta
    gb_seq = SeqIO.read(gb_file, "genbank")
    # print(gb_seq.seq)
    complete_seq = str(gb_seq.seq)
    # print(gb_seq)
    # print(gb_seq.annotations["accessions"])

    complete_ana = ">" + gb_seq.id + ":" + \
        gb_seq.annotations["accessions"][0] + " " + gb_seq.description + "\n"

    complete_fasta = format_fasta(complete_ana, complete_seq, 70)

    # print(type(gb_seq.features))

    # 提取 CDS 序列并格式为 fasta
    cds_fasta = ""

    for ele in gb_seq.features:
        if ele.type == "CDS":  # or ele.type == "gene":
            # print(ele.qualifiers)
            cds_seq = ""
            tmp_list = []
            # print(ele.location.parts)
            for ele1 in ele.location.parts:
                # print(ele1.start)
                # print(ele1.end)
                #print(int(re.findall(r'\d+', str(ele1.start))[0]))
                tmp_list.append(int(re.findall(r'\d+', str(ele1.start))[0]))
                tmp_list.append(int(re.findall(r'\d+', str(ele1.end))[0]))
                cds_seq += complete_seq[ele1.start:ele1.end]
            cds_ana = ">" + gb_seq.id + \
                " [" + str(tmp_list[0]+1)+".." + str(tmp_list[-1])+"]" + \
                " [gene=" + ele.qualifiers['gene'][0] + "]" + "\n"
            print(cds_ana)
            cds_fasta += format_fasta(cds_ana, cds_seq, 70)
            if (f_cds):
                break
    return cds_fasta, complete_fasta


if __name__ == '__main__':
    # 文件输出路径
    cds_file = "out/cds.fasta"
    complete_file = "out/complete.fasta"
    # genbank 文件路径
    res_dir = "F:\\Hibiscus_sabdariffa\\111"  # 改 改 改
    cds_file_obj = open(cds_file, "w")
    complete_file_obj = open(complete_file, "w")
    for file in os.listdir(res_dir):
        # print(os.sep)
        cds_fasta, complete_fasta = get_cds(res_dir + os.sep + file, False)
        cds_file_obj.write(cds_fasta)
        complete_file_obj.write(complete_fasta)
