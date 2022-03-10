#!/usr/bin/python3
# -*- coding : utf-8 -*-
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   根据第5套密码子表翻译DNA')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input', type=str, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

coding_dna = Seq(args.input)
print(coding_dna.translate(table=5))
