#!/usr/bin/env python3
import sys
import argparse

def main():
    arg_parser = argparse.ArgumentParser( description = 'Transform output friom IUPred and moreRONN into a format that is more easily-parseable.' )

    arg_parser.add_argument( '--input', '-i', help = "The input file to transform. This should be a file either output from IUPred or moreRONN. "
                                                     "It will be inferred from which program the output is. ",
                             type = str
                           )
    arg_parser.add_argument( '--output', '-o', help = "Name of file to write output to. File will contain a header." )


    args = arg_parser.parse_args() 

if __name__ == '__main__':
    main()
