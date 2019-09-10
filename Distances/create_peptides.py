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
    argp.add_argument( '--output', help = "Name of fasta file to write random sequences to.",
                      default = 'random.fasta'
                    )

    args = argp.parse_args()

    with open( args.output, 'w' ) as out_file:
        for index in range( args.num_seqs ):
            seq = generate_nt_seq( args.size )

            out_file.write( f'>seq_{index}\n{seq}\n' )

def generate_nt_seq( size ):
    seq = ''

    # codon size will not change
    codon_size = 3

    for index in range( size // codon_size ):
        seq += generate_non_stop_codon()
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
