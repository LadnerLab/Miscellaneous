#!/usr/bin/env python

import optparse

#This script takes as input multiple, multiple sequence alignments in fasta format, 
#     concatenates them all together and creates the necessary file to use --merge with MAFFT

def main():
    usage = '%prog [options]  msa1  [msa2 ...]'
    p = optparse.OptionParser()
    p.add_option('-o', '--out',  help='Base string for output files. [None, REQ]')
    opts, args = p.parse_args()
    
    totalN=[]
    totalS=[]
    
    with open("%s_msaKey.txt" % (opts.out),  "w") as fout:
        count=1
        for each in args:
            n,s = read_fasta_lists(each)
            fout.write("%s  # %s\n" % (" ".join([str(x) for x in range(count,count+len(n))]), each))
            totalN+=n
            totalS+=s
            count+=len(n)
    write_fasta(totalN, totalS, "%s.fasta" % (opts.out))
        
#----------------------End of main()

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
