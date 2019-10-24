#!/usr/bin/env python

import optparse, re

#Pull out a subset of sequences based on  taxonomy IDs included in names

def main():
    usage = '%prog [options]'
    p = optparse.OptionParser()
    p.add_option('-f', '--fasta',  help='Aligned fasta. [None, REQ]')
    p.add_option('-o', '--out',  help='out name. [None, REQ]')
    p.add_option('-i', '--ids', help='Taxonomic IDs of interest, comma-separated [None, REQ]')
    p.add_option('-t', '--tax', help='Taxonomic categories for IDS, comma-separated in same order as --ids. Options=fm,gn,sp. [None, REQ]')
    opts, args = p.parse_args()
    
    idD = parseTax(opts.ids.split(","), opts.tax.split(","))
    subfasta(idD, opts)
        
#----------------------End of main()

def parseTax(ids, tax):
    idD = {}
    code={"sp":0, "gn":1, "fm":2}
    for i,each in enumerate(ids):
        idD[each] = code[tax[i]]
    return(idD)
        

def taxMatch(name, idD):
    taxInfo = parse_tax(name)
    if taxInfo:
        for id, type in idD.items():
            if taxInfo[type] == id:
                return True
        return False
    else:
        return False
    
def parse_tax(name):
    oxpat = re.compile("OXX=(\d*),(\d*),(\d*),(\d*)")
    tax_id = oxpat.search(name)
    if tax_id:
        species,genus,family = (tax_id.group(2),tax_id.group(3),tax_id.group(4))
        return species,genus,family
    else:
        #print(name)
        return None

def subfasta(idD, opts):
    #Read in seqs
    names, seqs = read_fasta_lists(opts.fasta)
    
    new_names=[] 
    new_seqs=[]
    #Step through each sequence, check the amount of missing data and decide whether to keep it
    for i, s in enumerate(seqs):
        if taxMatch(names[i],idD):
            new_names.append(names[i])
            new_seqs.append(s)

    if len(new_names)>0: write_fasta(new_names, new_seqs, opts.out)
    else: print ('!!!!! No seqs met criteria !!!!!')
    
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

