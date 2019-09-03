#!/usr/bin/env python3
import argparse

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

    # parse sequences from fasta, grab IDs for each


    # find those that have changed (ID remains, OX different)

    # filter out IDs that are not from the fasta (if included)

    # write the output

def parse_records( fname ):
    records = list()

    new_record = lambda x: x[ 0 ] == 'ID'
    record_id  = lambda x: x[ 0 ] == 'OX'
    get_id     = lambda x: x.split( '=' )[ 1 ][ :-1 ]

    with open( fname, 'r' ) as of:
        for line in of:
            spl = line.strip().split()

            try:
                
                if new_record( spl ):
                    records.append( DataRecord( id = spl[ 1 ] ) )
                elif record_id( spl ):
                    records[ -1 ].set_ox( get_id( spl[ 1 ] ) )
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
