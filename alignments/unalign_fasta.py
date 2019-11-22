#!/usr/bin/env python

# Written by Jason Ladner

from __future__ import division
import optparse, re, os, itertools

#This script removes all '-' characters from each sequence. 
#The idea is to take an aligned fasta and convert to an unaligned fasta

def main():
    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
#Inputs and Outputs
    p.add_option('-f', '--fasta',  help='Aligned fasta with ambiguous characters. [None]')

    opts, samples = p.parse_args()
    
    names, seqs = read_fasta_lists(opts.fasta)
    new_seqs = [x.replace('-', "") for x in seqs]
    write_fasta(names, new_seqs, "unalign_%s" % (opts.fasta))

#End of main------------------------------------------


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

###-------------->>>

if __name__ == '__main__':
    main()    



