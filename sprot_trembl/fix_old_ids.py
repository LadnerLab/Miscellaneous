#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo
import re
import pandas as pd

def main():
    argp = argparse.ArgumentParser( description = 'Update old ids in a name map with the new ones.'
                                  )

    argp.add_argument( '--fixed', help = 'Name of the file that corrects ids (output by find_changes.py)' )
    argp.add_argument( '--original', help = "Name of the map file to fix. This file should be a tab-delimited "
                       "file with two columns. The first column should be the original name, the second the "
                       "name of the encoded peptide."
                     )
    argp.add_argument( '--output', help = "Name of file to output. Format of the output "
                       "is the same as that of the input, but the names will be fixed to include "
                       "the correct species-level IDs.", default = "output.tsv"
                     )

    args = argp.parse_args()

    id_map = get_new_ids( args.fixed )

    with open( args.original, 'r' ) as orig_map, open( args.output, 'w' ) as new_map:
        for line in orig_map:
            spl = line.strip().split( '\t' )

            id = get_id( spl[ 0 ] )
            out_id = spl[ 0 ]

            if id in id_map:
                out_id = replace_id( spl[ 0 ], id_map[ id ] )
            new_map.write( f'{out_id}\t{spl[ 1 ]}\n' )

def get_new_ids( fname ):
    output = dict()

    with open( fname, 'r' ) as of:
        for lineno, line in enumerate( of ):
           if lineno: # header skip
               spl = line.strip().split( '\t' )
               output[ spl[ 0 ] ] = spl[ 2 ]
    return output

def get_id( name ):
    id_re = re.compile( "ID=([\S]+)\s[\S]+\sOXX=\d*,\d*,\d*,\d*" )

    search = re.search( id_re, name )
    if search:
        return search.group( 1 )
    return None

def replace_id( old_name, new_id ):
    id_re = re.compile( "(ID=[\S]+\s[\S]+\sOXX=)(\d*),(\d*),(\d*),(\d*)" )

    return re.sub( id_re, '\\1' + '\\2' + ',' + new_id + ',' + '\\4' + ',' + '\\5', old_name )
            

if __name__ == '__main__':
    main()
