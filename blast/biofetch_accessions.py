#!/usr/bin/env python

#from Bio import Entrez, SeqIO
from Bio import Entrez
import optparse

def main():

    #To parse command line
    usage = "usage: %prog [options]"
    p = optparse.OptionParser(usage)
    
    p.add_option('-e', '--email', default="jason.ladner@nau.edu", help='your email address [jason.ladner@nau.edu]')
    p.add_option('-d', '--db', default="nuccore", help='name of genbank database [nuccore]')
    p.add_option('-i', '--info', help='Tab-delimited file with info on seqs to obtain and names for output [None, REQ]')
#    p.add_option('-o', '--out', help='Name for output file [None, REQ]')
#    p.add_option('-b', '--batchsize', default=100, help='batch size for queries. [100]')
#    p.add_option('-r', '--retmax', default=10**9, help='retmax. [10**9]')

    opts, args = p.parse_args()

    Entrez.email=opts.email
    
    with open(opts.info, "r") as fin:
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            outName = "%s_%s_%s.faa" % (cols[0], cols[1], cols[2])
            with open(outName, "w") as fout:
                thishandle = Entrez.efetch(db=opts.db, id=cols[2], rettype="fasta_cds_aa", retmode="text")
                fout.write(thishandle.read())

#------------------------------------------------------------------>


       
###---------------------->>>

if __name__ == '__main__':
    main()
    

