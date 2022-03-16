from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio import SearchIO
<<<<<<< Updated upstream
xml = SearchIO.parse('F:/ref_tre/gene/blast/blastn.result.xml', 'blast-xml')
SearchIO.write(xml, 'F:/ref_tre/gene/blast', 'blast-tab')


# 只提供序列意味着BLAST会自动分配给你一个ID。你可能更喜欢用 SeqRecord 对象的format方法来包装一个fasta字符串，因为这个对象会包含fasta文件中已有的ID

record = SeqIO.read("m_cold.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))
=======
from Bio.Blast.Applications import NcbiblastxCommandline
blastx_cline = NcbiblastxCommandline(
    query="opuntia.fasta", db="nr", evalue=0.001, outfmt=5, out="opuntia.xml")
print(blastx_cline)
# (stdout, stderr) = blastx_cline()11111

"""
xml = SearchIO.parse('F:/ref_tre/gene/blast/blastn.result.xml', 'blast-xml')
SearchIO.write(xml, 'F:/ref_tre/gene/blast', 'blast-tab')

#只提供序列意味着BLAST会自动分配给你一个ID。你可能更喜欢用 SeqRecord 对象的format方法来包装一个fasta字符串，因为这个对象会包含fasta文件中已有的ID
record = SeqIO.read("m_cold.fasta", format="fasta")
result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"))


#这里我们需要注意下：因为用 result_handle.read() 来读取BLAST结果只能用一次 - 再次调用 result_handle.read() 会返回一个空的字符串.
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()

#这些做好后，结果已经存储在 my_blast.xml 文件中了并且原先的handle中的数据 已经被全部提取出来了(所以我们把它关闭了)。但是，BLAST解析器的 parse 函数（描述见 7.3) 采用一个文件句柄类的对象，所以我们只需打开已经保存的文件作为输入。
result_handle=open("my_blast.xml")
"""
>>>>>>> Stashed changes
