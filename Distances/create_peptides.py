#!/usr/bin/env python3
import argparse
import random

def main():
    argp = argparse.ArgumentParser( description = "Create some number of randomly-generated nt kmers" )

    argp.add_argument( '--num_seqs', help = "The number of sequences to generate.",
                       type = int, default = 244000
                     )
    argp.add_argument( '--size', help = "The size of sequence to generate.",
                       type = int, default = 90
                     )
    argp.add_argument( '--codons', help = "Name of file to generate codons from." )
    argp.add_argument( '--output', help = "Name of fasta file to write random sequences to.",
                      default = 'random.fasta'
                    )

    args = argp.parse_args()

    if args.codons:
        codons = get_codons( args.codons )

    with open( args.output, 'w' ) as out_file:
        for index in range( args.num_seqs ):
            if args.codons:
                seq = generate_nt_seq( args.size, lambda: random.choice( codons ) )
            else:
                seq = generate_nt_seq( args.size, generate_non_stop_codon )

            out_file.write( f'>seq_{index}\n{seq}\n' )

def get_codons( fname ):
    out = list()
    with open( fname, 'r' ) as of:
        for line in of:
            spl = line.strip().split( ',' )
            out.append( spl[ 1 ] )
    return out

def generate_nt_seq( size, gen_fn ):
    seq = ''

    # codon size will not change
    codon_size = 3

    for index in range( size // codon_size ):
        seq += gen_fn()
    return seq
    
def generate_non_stop_codon():

    stop_codons = { 'TAG', 'TAA', 'TGA' }
    cod = generate_codon()
    while cod in stop_codons:
        cod = generate_codon()
    return cod

def generate_codon():
    nt_map = { 1: 'A', 2: 'G', 3:'T', 4:'C' }
    return ''.join( map( nt_map.get, gen_n_random_values( 3, 1, 4 ) ) )

def gen_n_random_values( n, a, b ):
    out = list()

    for index in range( n ):
       out.append( random.randint( a, b ) )
    return out


if __name__ == '__main__':
    main()
