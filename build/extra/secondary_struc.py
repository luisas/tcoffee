#!/usr/bin/env python

import os
import sys
import re
import string


seq_file = sys.argv[1]
pdb_file = sys.argv[2]
out_file = sys.argv[3]

chain = pdb_file[-5]

regular_ex = re.compile(r"ATOM\s*\S+\s+\S+\s+(\S+)\s+([A-Za-z])\s*(\d+)\s+")
regular_ex2 = re.compile(r"^[^:]+:(\d+)-[^:]+:(\d+)")



def read_fasta(file_name):
	
	f_file = file(file_name)
	h_seq = {}
	while True:
		line = f_file.readline()
		if len(line) == 0:
			break
		if line[-1] == '\n':
			line = line[:-1]
		if (len(line)>0):
			if (line[0] == ">"):
				name = line[1:]
				h_seq[name] = ""
			else:
				h_seq[name] += (line)
	f_file.close()
	return h_seq



#make hash from 
def pdb2seq(pdb_name, chain):
	
	pdb_f = open(pdb_name)
	seq_num = 0
	old_res = -1
	pdb_h = {}
	seq = []
	chain_found = False
	while True:
		line = pdb_f.readline()
		m = line[:4]
		
		if m == "ATOM":
			
			
			if line[21] == chain:
				chain_found = True
				m3 = line[22:26].strip()
				if (m3 != old_res):
					seq.append(line[19])
					pdb_h[m3] = seq_num
					old_res = m3
					seq_num += 1
		if chain_found & ((len(line)==0) | (line[:3] == "TER")):
			break
	pdb_f.close()
	
	return (pdb_h, string.join( seq, '' ))



def make_hash_seq2seq(aln_seq1, aln_seq2):
	h_match = {}
	length = len(aln_seq1)
	pos_seq1 = 0;
	pos_seq2 = 0;
	for aln_pos in range(0, length):
		if (aln_seq1[aln_pos] != '-') & (aln_seq2[aln_pos] != '-'):
			h_match[pos_seq1] = pos_seq2
		if (aln_seq1[aln_pos] != '-'):
			pos_seq1 += 1
		if (aln_seq2[aln_pos] != '-'):
			pos_seq2 += 1
	return h_match



def ungap_seq(seq):
	reg = re.compile(r"-")
	s = reg.sub("", seq)
	return s






seq_a = []
(pdb_h,seq) = pdb2seq(pdb_file,chain)






seq_f = open(seq_file,"a+")
seq_name = seq_f.readline()[1:-1]
seq_f.write(">PDB\n"+seq+"\n")
seq_f.close()



os.system("t_coffee -method fast_pair -in "+ seq_file + " -outfile " + seq_name +"out -output fasta 2>/dev/null >/dev/null")


align_h = read_fasta(seq_name+"out")



pdb_to_seq_h = make_hash_seq2seq(align_h["PDB"], align_h[seq_name])




home_path  =  os.environ.get('HOME') + '/'


os.system("rnass.py " + pdb_file + " " + chain + " >/dev/null 2>/dev/null")



sec_file = open(pdb_file+"."+chain+'.pairs')
sec_file_out = open(out_file,'w')



sec_file_out.write("! TC_LIB_FORMAT_01\n")
sec_file_out.write("\n")
sec_file_out.write("1\n")
sec_file_out.write(seq_name + " " + str(len(align_h[seq_name])) + " " + ungap_seq(align_h[seq_name]) + "\n")
sec_file_out.write("!generated by secondary_struc.py at some day and some day\n")
sec_file_out.write("#1 1\n")



while True:
	line = sec_file.readline()
	if len(line) == 0:
		break
	m = regular_ex2.match(line)
	if m != None:
		m1 = pdb_h[m.group(1)]
		m2 = pdb_h[m.group(2)]
		if (m1 in pdb_to_seq_h) & (m2 in pdb_to_seq_h):
			sec_file_out.write(str(pdb_to_seq_h[m1]+1) + " " + str(pdb_to_seq_h[m2]+1) + " 100\n")

		
sec_file_out.write("! SEQ_1_TO_N\n")
sec_file.close()
sec_file_out.close()

os.remove(pdb_file+"."+chain+'.pairs')
os.remove(pdb_file+"."+chain+'.ss')

