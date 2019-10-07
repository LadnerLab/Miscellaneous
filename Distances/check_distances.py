#!/usr/bin/env python3
import sys
import protein_oligo_library as oligo

def main():
    in_dist  = sys.argv[ 1 ]

    with open( in_dist, 'r' ) as of:
        for lineno, line in enumerate( of ):
            if lineno:
                s1, dist, s2 = line.strip().split( '\t' )
                dist = int( dist )
                assert( distance( s1, s2 ) == dist ) 
    print( "All checks passed!" )


def distance( a, b ):
    return sum( [ a_c != b_c for a_c, b_c in zip( a, b ) ] )

if __name__ == '__main__':
    main()

