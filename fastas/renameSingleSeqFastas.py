#!/usr/bin/env python

import optparse

#Rename seq in fasta to match name of file

def main():
    usage = '%prog [options] fasta1 [fasta2 ...]'
    p = optparse.OptionParser()
    opts, args = p.parse_args()
    
    for f in args:
        n,s = read_fasta_lists(f)
        if len(n)>1:
            print("This script is only meant to be used for fastas containing a single sequence each. Skipped: %s" % (f))
        else:
            write_fasta([".".join(f.split(".")[:-1])], s, f)
        
#----------------------End of main()

    
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

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

