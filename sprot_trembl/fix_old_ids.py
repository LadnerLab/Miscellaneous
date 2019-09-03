#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo
import re

def main():
    argp = argparse.ArgumentParser( description = 'Update old ids in a name map with the new ones.'
                                  )

    argp.add_argument( '--fixed', help = 'Name of the file that corrects ids (output by find_changes.py)' )
    argp.add_argument( '--original', help = "Name of the file to fix." )
    argp.add_argument( '--output', help = "Name of file to output.", default = "output.tsv" )

    args = argp.parse_args()



if __name__ == '__main__':
    main()
