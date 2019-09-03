#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo
import re

def main():
    argp = argparse.ArgumentParser( description = 'Find changes in Taxonomic ID from one version of '
                                    'sprot/tremble data files to another. Optionally '
                                    'only report results whose ids are in the sequences in '
                                    'a fasta'
                                  )

    argp.add_argument( '--first', help = "The first data file to check." )
    argp.add_argument( '--second', help = "The second data file to check." )
    argp.add_argument( '--filter', help = "Fasta file of sequences to include. " )
    argp.add_argument( '--output', help = "Name of file to output." )

    args = argp.parse_args()

    # parse sequences from the first and second (only need ID and OX)
    first_recs  = parse_records( args.first )
    second_recs = parse_records( args.second )

    # slightly trim search space
    to_search = first_recs.keys() & second_recs.keys()

    # parse sequences from fasta, grab IDs for each
    fasta_records = get_records_from_fasta( args.filter )

    # find those that have changed (ID remains, OX different)
    changed_ids = filter( lambda x: first_recs[ x ] != second_recs[ x ],
                          to_search
                        )

    # filter out IDs that are not from the fasta (if included)
    if fasta_records:
        target_ids = set( changed_ids ) & fasta_records.keys() 
    else:
        target_ids = set( changed_ids )

    # write the output

def get_records_from_fasta( fname ):
    records = dict()
    line_re = re.compile( "ID=([\S]+)\s[\S]+\sOXX=\d*,(\d*),\d*,\d*" )
    if fname:
        names, seqs = oligo.read_fasta_lists( fname )
        for name in names:
            name_search = re.search( line_re, name )
            if name_search:
                id = name_search.group( 1 )
                ox = name_search.group( 2 )

                # records.append( DataRecord( id = id, ox = ox ) )
                records[ id ] = str( ox )
            else:
                print( "No match: ", name )
    return records

def parse_records( fname ):
    records = dict()

    new_record = lambda x: x[ 0 ] == 'ID'
    record_id  = lambda x: x[ 0 ] == 'OX'
    get_id     = lambda x: x.split( '=' )[ 1 ][ :-1 ]

    with open( fname, 'r' ) as of:
        for line in of:
            spl = line.strip().split()
            try:
                
                if new_record( spl ):
                    # records.append( DataRecord( id = spl[ 1 ] ) )
                    records[ spl[ 1 ] ] = ''
                    last_id = spl[ 1 ]
                elif record_id( spl ):
                    records[ last_id ] = get_id( spl[ 1 ] )
            except IndexError:
                print( line )
                
    return records



class DataRecord:
    def __init__( self, id = '', ox = '' ):
        self.id = id
        self.ox = ox

    def get_id( self ):
        return self.id

    def get_ox( self ):
        return self.ox

    def set_id( self, new_id ):
        self.id = new_id

    def set_ox( self, new_ox ):
        self.ox = new_ox

    def __hash__( self ):
        return hash( self.id )

    def __eq__( self, other ):
        return self.id == other.id 

    def __ne__( self, other ):
        return not self.__eq__( other )

if __name__ == '__main__':
    main()
