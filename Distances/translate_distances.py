#!/usr/bin/env python3
import argparse
from enum import Enum

class CountMode( Enum ):
    NT = 0
    AA = 1


def main():
    argp = argparse.ArgumentParser( description = "Given the nearest neighbors of either NT/AA sequences and a map "
                                                  "associating AA <-> NT sequences, determine the translated "
                                                  "distances."
                                  )
    argp.add_argument( '-i', '--input', help = "Input distances to translate. Seqeuence type (NT or AA ) will be inferred." )
    argp.add_argument( '-m', '--map', help = "Map associating NT and AA sequences." )
    argp.add_argument( '-o', '--output', help = "Name of the translated file to write output to." ) 

    args = argp.parse_args()

    input_distances = parse_distances( args.input )
    input_format = nt_or_aa( input_distances )
    input_map = parse_map( args.map, fmt = input_format )

    output_distances = create_distances( input_distances,
                                         input_map
                                       )

def create_distances( in_dist, in_map ):
    output = list()
    hamm = lambda x, y: sum( [ a != b for a, b in zip( x , y ) ] )

    for distance in in_dist:
        seq, dist, neighbor = distance

        tran_seq = in_map[ seq ]
        tran_neighbor = in_map[ neighbor ]

        dist = hamm( tran_seq, tran_neighbor )
        new = tuple( tran_seq, str( dist ), tran_neighbor )

        output.append( new )
    return output

def parse_map( filename, fmt = CountMode.AA ):
    # input is formatted: NT\tAA, so this allows us to select
    # the necessary index
    key_index = fmt == CountMode.AA
    value_index = not( key_index )
    output = dict()
    
    with open( filename, 'r' ) as open_f:
        for line in open_f: # no header
            str_spl = line.strip().split( '\t' )

            new_key = str_spl[ key_index ]
            new_val = str_spl[ val_index ]

            output[ new_key ] = new_val 
    return output

def parse_distances( ifname ):
    output = list()

    with open( ifname, 'r' ) as open_f:
        for lineno, line in enumerate( ifname ):
            if lineno: 
                spl_str = line.strip().split( '\t' )

                output.append( tuple( spl_str ) )
                
    return output

def nt_or_aa( input_tuples, num_check = 10 ):
    aa_count = 0
    nt_count = 0

    aa = lambda x: set( x ) == set( [ 'A', 'T', 'G', 'C' ] )
    
    for x in range( num_check ):
        current_tup = input_tuples[ x ]

        if aa( current_tup[ 0 ] ):
            aa_count += 1
        else:
            nt_count += 1

        if aa( current_tup[ 1 ] ):
            aa_count += 1
        else:
            nt_count += 1

    if aa_count > nt_count:
        return CountMode.AA

    return CountMode.NT




if __name__ == '__main__':
    main()
