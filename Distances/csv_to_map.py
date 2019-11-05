#!/usr/bin/env python3
import sys

def main():

    if len( sys.argv ) != 3:
        print( f"USAGE: { sys.argv[0] } input_csv output_tsv" )
        sys.exit( 1 )

    in_csv  = sys.argv[ 1 ] 
    out_tsv = sys.argv[ 2 ]

    with open( in_csv, 'r' ) as in_f, open( out_tsv, 'w' ) as out_f:
        for lineno, line in enumerate( in_f ):
            if lineno:
                spl = line.strip().split( ',' )
                try:
                    peptide    = spl[ 1 ]
                    nucleotide = spl[ 2 ]

                    out_f.write( f'{ nucleotide }\t{ peptide }\n' )
                except IndexError:
                    print( f"Error on line {lineno}: {line}" )


if __name__ == '__main__':
    main()
