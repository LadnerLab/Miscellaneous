#!/usr/bin/env python3
import argparse
import fastatools

def main():
    argparser = argparse.ArgumentParser( "Output a FASTA containing only unique sequences." )
    argparser.add_argument( '-i','--input', help = "Name of input file." )
    argparser.add_argument( '-o','--output', help = "Name of output file This file will contain the "
                                                    "same sequences as the input file, but duplicates will not be included."
                         )

    args = argparser.parse_args()

    in_names, in_seqs = fastatools.read_fasta_lists( args.input )

    out_names, out_seqs = list(), list()
    seen_seqs = set()

    for name, seq in zip( in_names, in_seqs ):
        if seq not in seen_seqs:
            seen_seqs.add( seq )
            out_names.append( name )
            out_seqs.append( seq )

    fastatools.write_fasta( out_names, out_seqs, args.output )

if __name__ == '__main__':
    main()
