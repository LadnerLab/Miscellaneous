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

    # parse sequences from fasta, grab IDs for each


    # find those that have changed (ID remains, OX different)

    # filter out IDs that are not from the fasta (if included)

    # write the output

class DataRecord:
    def __init__( self, id = '', ox = '' ):
        self.id = id
        self.ox = ox

    def __hash__( self ):
        return hash( self.id )

    def __eq__( self, other ):
        return self.id == other.id 

    def __ne__( self, other ):
        return not self.__eq__( other )

if __name__ == '__main__':
    main()
