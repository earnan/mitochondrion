#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   mt_translation_dna.py
#         Author:   yujie
#    Description:   mt_translation_dna.py
#        Version:   1.0
#           Time:   2022/03/14 10:07:30
#  Last Modified:   2022/03/14 10:07:30
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import os
import re

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   根据第n套密码子表翻译DNA')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input', type=str, required=False)
optional.add_argument('-n', '--number', type=int, default=5, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

coding_dna = Seq(args.input)
print(coding_dna.translate(table=args.number))
