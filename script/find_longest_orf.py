import os
import argparse
from Bio import SeqIO
import re

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for finding longest orf and its amino acid sequence from sequence in fasta')

parser.add_argument('--input', '-i', metavar='seq.fasta',
        help='Input fasta',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='result.fasta',
        help='Output file name where resultant amino acid sequence will be written to')

args = parser.parse_args()

record_iterator = SeqIO.parse(args.input, "fasta")

record = next(record_iterator)
id = record.id
seq = str(record.seq)

orf = max(re.findall(r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)',seq), key = len)

orf_obj = Seq(orf)
aa = orf_obj.translate()

SeqIO.write(aa, args.output, "fasta")

#print(str(aa.seq))