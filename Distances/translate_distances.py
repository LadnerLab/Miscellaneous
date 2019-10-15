#!/usr/bin/env python3
import argparse

def main():
    argp = argparse.ArgumentParser( description = "Given the nearest neighbors of either NT/AA sequences and a map "
                                                  "associating AA <-> NT sequences, determine the translated "
                                                  "distances."
                                  )
    argp.add_argument( '-i', '--input', help = "Input distances to translate. Seqeuence type (NT or AA ) will be inferred." )
    argp.add_argument( '-m', '--map', help = "Map associating NT and AA sequences." )
    argp.add_argument( '-o', '--output', help = "Name of the translated file to write output to." ) 

    args = argp.parse_args()

if __name__ == '__main__':
    main()
