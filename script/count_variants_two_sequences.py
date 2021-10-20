import os
import argparse
from Bio import SeqIO

def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    return arg

parser = argparse.ArgumentParser(description= 'Script for count the number of variants between two sequences from MSA output')

parser.add_argument('--input', '-i', metavar='seqs.afa',
        help='Input MSA output',
        type=lambda x: is_valid_file(parser, x))

parser.add_argument('--output', '-o', metavar='result.txt',
        help='Onput file name where result will be written to')

parser.add_argument('--header', dest='header', default=False, action='store_true')

args = parser.parse_args()

record_iterator = SeqIO.parse(args.input, "fasta")

#First record will be regarded as reference
first_record = next(record_iterator)
first_id = first_record.id
first_seq = str(first_record.seq)

second_record = next(record_iterator)
second_id = second_record.id
second_seq = str(second_record.seq)

assert len(first_seq) == len(second_seq)
aln_len = len(first_seq)

count_dict = { "match" : 0,
               "mismatch" : 0,
               "insertion" : 0,
               "deletion" : 0,
               "low_cov" : 0}
indel_ongoing = False

for i in range(len(first_seq)):
    first_seq_chr = first_seq[i]
    second_seq_chr = second_seq[i]
    
    if first_seq_chr == "N" or second_seq_chr == "N": # Hard masked region should not be considered
        count_dict["low_cov"] += 1
        continue
    
    if first_seq_chr == second_seq_chr:
        count_dict["match"] += 1
    
    else:
                    
        if indel_ongoing:
            if first_seq_chr != "-" and second_seq_chr != "-":
                count_dict["mismatch"] += 1
                indel_ongoing = False
                
            else:
                continue
                
        else:
            if first_seq_chr != "-" and second_seq_chr != "-":
                count_dict["mismatch"] += 1

            else:
                if first_seq_chr == "-":
                    count_dict["insertion"] += 1
                    indel_ongoing = True
                else:
                    count_dict["deletion"] += 1
                    indel_ongoing = True

outfile = open(args.output, "w")
header_list = ["sample", "aln_len", "match", "mismatch", "insertion", "deletion", "low_cov"]
header = '\t'.join(header_list) + "\n"

base = os.path.basename(args.input)
sample_name = os.path.splitext(base)[0]
result_list = [sample_name, str(aln_len), str(count_dict["match"]), str(count_dict["mismatch"]), str(count_dict["insertion"]), str(count_dict["deletion"]), str(count_dict["low_cov"])]
result = '\t'.join(result_list) + "\n"

if args.header:
    outfile.write(header)
outfile.write(result)
outfile.close()

