#!/usr/bin/env python3
import sys
import argparse
import re

def main():
    arg_parser = argparse.ArgumentParser( description = 'Transform output friom IUPred and moreRONN into a format that is more easily-parseable.' )

    arg_parser.add_argument( '--input', '-i', help = "The input file to transform. This should be a file either output from IUPred or moreRONN. "
                                                     "It will be inferred from which program the output is. ",
                             type = str
                           )
    arg_parser.add_argument( '--output', '-o', help = "Name of file to write output to. File will contain a header." )


    args = arg_parser.parse_args() 

    parser, transform_fn = get_parser( args.input )

    transformed_out = list()

    with open( args.input, 'r' ) as of:
        for record in parser.parse( of ):
            transformed_rec = transform_fn( record )
            transformed_out.append( transformed_rec )
        


def get_parser( input_fname ):
    parser, transformer = MoronParser(), MoronTransformer().get_transform()

    if not parser.peek( input_fname ):
        parser, transformer = IUPredParser(), IUPredTransformer().get_transform()

        if not parser.peek( input_fname ):
            raise ValueError( "Input file is not in either IUPred or MoreRONN format" )  

    return parser, transformer

class MoronRecord:
    def __init__( self, seq_id = "", seq = "", score = "", pos_scores = None ):
        if not pos_scores:
            self.pos_scores = list()
        else:
            self.pos_scores = pos_scores

        self.seq_id = seq_id
        self.seq = seq
        self.score = score

class MoronTransformer:
    def transform( self, moron_rec ):
        return self._transform( moron_rec )

    def _transform( self, record ):
        output = '%s\t' % ( record.seq_id )

        scores = list()

        for aa, score in record.pos_scores:
            scores.append( score )

        return output + ','.join( scores )

    def get_transform( self ):
        return self.transform

class IUPredTransformer:
    def transform( self, str_list ):
        return _transform( str_list )

    def get_transform( self ):
        return self.transform

class MoronParser:

    def parse( self, input_file ):
        record = None
        last_line = ""

        is_identifier = lambda x: '(Sequence ' in x
        is_aa_score = lambda x: len( x.split( '\t' ) ) == 2
        
        for line in input_file:
            if is_identifier( line ):
                if record:
                    yield record
                record = MoronRecord( seq_id = line.split( ':' )[ 1 ].strip() )

            elif is_aa_score( line ):
                split = line.strip().split( '\t' )
                record.pos_scores.append( tuple( split ) )

            last_line = line

    def peek( self, fname ):
        first_re = re.compile( '\(Sequence (\d+) of (\d+)\):' )

        with open( fname, 'r' ) as of:
            for line in of:

                if len( line ) and re.match( first_re, line.strip() ):
                    return True
            return False


class IUPredParser:
    def parse( self, input_file ):
        pass

    def peek( self, fname ):
        with open( fname, 'r' ) as of:
            for line in of:
                if line[ 0 ] == '>':
                    return True
                return False
        return False

if __name__ == '__main__':
    main()
