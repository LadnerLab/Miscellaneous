#!/usr/bin/env python3
import argparse

def main():
    args = argparse.ArgumentParser( description = 'Find changes in Taxonomic ID from one version of '
                                    'sprot/tremble data files to another. Optionally '
                                    'only report results whose ids are in the sequences in '
                                    'a fasta'
                                  )

    args.add_argument( '--first', help = "The first data file to check." )
    args.add_argument( '--second', help = "The second data file to check." )
    args.add_argument( '--filter', help = "Fasta file of sequences to include. " )
    args.add_argument( '--output', help = "Name of file to output." )

    # parse sequences from the first and second (only need ID and OX)

    # parse sequences from fasta, grab IDs for each


    # find those that have changed (ID remains, OX different)

    # filter out IDs that are not from the fasta (if included)

    # write the output

if __name__ == '__main__':
    main()
