#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo

def main():
    argp = argparse.ArgumentParser( description = "Given a map of codons to Amino Acids translate a "
            "nt fasta file to an aa fasta file." 
                                  )
    argp.add_argument( '--fasta', help = "Input fasta to parse." )
    argp.add_argument( '--map', help = "Name of file mapping codons to amino acids" )
    argp.add_argument( '--output', help = "Name of file to write faa file to.", default = "output.fasta" )

    args = argp.parse_args()

    orig_seqs = parse_fasta( args.fasta )
    aa_nt_map = parse_map( args.map )

    list( map( lambda x: translate( aa_nt_map, x ), orig_seqs  ) )

    write_output( args.output,  orig_seqs )

def write_output( fname, seqs ):
    with open( fname, 'w' ) as of:
       for s in seqs:
            of.write('>{s.name}\n{s.seq}\n')


def translate( aa_map, seq ):
    new_seq = ''
    for codon_idx in range( 0, len( seq.seq ), 3 ):
        try:
            cod = seq.seq[ codon_idx : codon_idx + 3 ]
            new_seq += aa_map[ cod ]
        except KeyError:
            print( seq.seq, cod, codon_idx )
    seq.seq = new_seq 

def parse_map( fname ):
    out = dict()

    with open( fname, 'r' ) as of:
        for line in of:
            spl = line.strip().split( ',' )

            out[ spl[ 1 ] ] = spl[ 0 ] 
    return out

def parse_fasta( fname ):
    n, s = oligo.read_fasta_lists( fname )

    return [ Sequence( na, se ) for na, se in zip( n, s ) ]

class Sequence:
    def __init__( self, name, seq ):
        self.name = name
        self.seq  = seq


if __name__ == '__main__':
    main()
