#!/usr/bin/env python

import optparse, os

# Checks the final length of one or more alignments
# Will also look for an unaligned version (based on naming assumption)
    # If it finds one, it will also report the difference betweeen the alignment length and the length of the longest unaligned seq

def main():
    usage = '%prog [options] align1.fasta [align2.fasta ...]'
    p = optparse.OptionParser(usage)
    opts, args = p.parse_args()
    
    for each in args:
        names, seqs = read_fasta_lists(each)
        lens = [len(s) for s in seqs]
        if len(set(lens))>1:
            print("%s is NOT properly aligned. Observed lengths:%s" % (each, [str(x) for x in set(lens)]))
        else:
            print(each)
            print("\tAlignment: %d" % (lens[0]))
            putUnaligned = "%s.%s" % ("_".join(each.split("_")[:-1]), each.split(".")[-1])
            if os.path.isfile(putUnaligned):
                uN, uS = read_fasta_lists(putUnaligned)
                uLens = [len(s) for s in uS]
                diff = lens[0]-max(uLens)
                print("\tVs Unaligned: %d (%.1f)" % (diff, diff/lens[0]*100))

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

###------------------------------------->>>>    

if __name__ == "__main__":
    main()

