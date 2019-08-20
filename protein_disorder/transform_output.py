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

    parser = get_parser( args.input )


def get_parser( input_fname ):
    parser = MoronParser()

    if moron_p.peek( input_fname ):
        pass
    else:
        parser = IUPredParser()

        if not parser.peek( input_fname ):
            raise ValueError( "Input file is not in either IUPred or MoreRONN format" )  

    return parser

class MoronParser:

    def parse( self, input_file ):
        pass

    def peek( self, fname ):
        first_re = re.compile( '\(Sequence (\d+) of (\d+)\):' )

        with open( fname, 'r' ) as of:
            for line in fname:

                if len( line ) and re.findall( first_re, line.strip() ):
                    return True
            return False


class IUPredParser:
    def parse( self, input_file ):
        pass

    def peek( self, fname ):
        pass

if __name__ == '__main__':
    main()
