#!/usr/bin/env python

import argparse, os
import fastatools as ft

from collections import defaultdict

#Pull out a subset of sequences based on  taxonomy IDs included in names

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--simple', action="store_true", default=False, help="Use this option if you want names in output fasta to be 'simple', not including full file path.")
    parser.add_argument('-n', '--numExp', type=int, help="Use this option if you want to limit the output to only sections that include all of your aligned genomes. If used, please provide the total number of aligned genomes.")

    reqArgs = parser.add_argument_group('Required Arguments')
    reqArgs.add_argument('-x', '--xmfa', help='XMFA input file', required=True)
    reqArgs.add_argument('-o', '--out',  help='Name for output fasta file.', required=True)

    args = parser.parse_args()
    
    fD = ft.read_fasta_dict_upper(args.xmfa)
    
    oD = defaultdict(str)
    subD = {}
    
    countsL = []
    
    for n,s in fD.items():
        nameParts = n.split()
        name = nameParts[-1]
        if args.simple:
            name = os.path.basename(name)
        
        if s[-1]=="=":
            s = s[:-1]
            subD[name]=s
            
            if args.numExp:
                if len(subD) == args.numExp:
                    for k,v in subD.items():
                        oD[k]+=v
                    countsL.append(len(subD))
            else:
                for k,v in subD.items():
                    oD[k]+=v
                countsL.append(len(subD))
            
            subD = {}
            
        else:
            subD[name]=s
    
    print(countsL)
        
    ft.write_fasta_dict(oD, args.out)
    
#----------------------End of main()


###------------------------------------->>>>    

if __name__ == "__main__":
    main()

