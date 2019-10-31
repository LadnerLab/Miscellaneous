#! /usr/bin/env python

from read_write_files import read_fasta_lists, write_fasta                 #Imports two functions, one that reads in name and sequences from a fasta file and the other writes out the new fasta in the correct format, both are defined in read_write_files.py, available at www.stanford.edu/~jtladner 
import sys

# This script extracts a subset of sequences from a phylip file and writes them to a new phylip file with the same name as the old file + a suffix that is specified on the command line
# The script accepts any number of phylip files for batch processing using the same subset of individuals

# Usage: python subset_fasta.py name_list prefix to_keep_or_remove file_1.fasta file_2.fasta...file_n.fasta
# Tip: '*.fasta' specifies all the fasta files in the current directory
# name_list should be a plain text file with a list of names specifying the subset of sequences you are interested in, each name should be on a new line
# Prefix can be any string, it will be appended to the new filenames before the '.phy' extension
# The variable 'to_keep_or_remove' should be either 0 (to remove the sequences specified in 'name_list') or 1 (to keep only the sequences in 'name_list')


def read_fasta_dict(file):
	names, seqs = read_fasta_lists(file)
	fasta_dict = dict(zip(names, seqs))
	return fasta_dict

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
	fin = open(file, 'r')
	count=0
	
	names=[]
	seqs=[]
	seq=''
	for line in fin:
		line=line.strip()
		if line and line[0] == '>':                #indicates the name of the sequence
			count+=1
			names.append(line[1:])
			if count>1:
				seqs.append(seq)
			seq=''
		else: seq +=line
	seqs.append(seq)
	
	return names, seqs

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
	fout=open(new_filename, 'w')
	for i in range(len(names)):
		fout.write(">%s\n%s\n" % (names[i], seqs[i]))
	fout.close()

if len(sys.argv) < 5: print ('Usage: subset_fasta.py name_list prefix to_keep_or_remove file_1.fasta file_2.fasta...file_n.fasta')


indiv = open(sys.argv[1], 'r')                  #Reads in the subset of individuals you are interested in
prefix=str(sys.argv[2])                         #Suffix to add onto the new filenames
to_keep_sub=int(sys.argv[3])                    #Should be either a 1 or 0. Specifies whether the subset in indiv variable should be kept or discarded
files=sys.argv[4:]                              #Phylips files
if to_keep_sub > 1 or to_keep_sub < 0:
	print ("Third argument should be either 1 or 0\n")


subset={}
#Creates a dictionary of all the names in name_list
for line in indiv:

	line = line.rstrip()
	cols = line.split('\t')
	subset[cols[0]]=""

#Steps through each phylip file and pulls out all the sequences that match the names in name_list
for file in files:
	IN = open(file, 'r')
	names, seqs=read_fasta_lists(file)
#	print names
	nametrim=[name.rstrip() for name in names]
	
	nopath_file=file.split('/')[-1]
	
	subnames=[]
	subseqs=[]
	
	for index in range(len(names)):
		if to_keep_sub:
			if nametrim[index] in subset:
				subnames.append(names[index])
				subseqs.append(seqs[index])

		else:
			if nametrim[index] not in subset:
				subnames.append(names[index])
				subseqs.append(seqs[index])
	write_fasta(subnames, subseqs, prefix + nopath_file)
	IN.close()

indiv.close()
