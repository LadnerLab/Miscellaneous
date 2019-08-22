#!/usr/bin/env python3
import sys
import argparse
import protein_oligo_library as oligo
import numpy as np


def main():
    argp = argparse.ArgumentParser( description = "Calculate mean scores (IUPred/MoreRONN) for a set of peptide probes, "
                                                  "report number of positions greater than a certain threshold."
                                  )

    argp.add_argument( '--probes', '-p', help = "FASTA file containing probes to check. "
                       "Positions from original sequence must be in name (generated by kmer_oligo)."
                     )
    argp.add_argument( '--iupred', help = "File containing transformed IUPred scores." )
    argp.add_argument( '--more_ronn', help = "File containing transformed MoreRONN scores." )
    argp.add_argument( '--output', '-o', help = "Name of file to write output to.", default = "probe_summaries.tsv" )
    argp.add_argument( '--thresholds', '-t', help = "Comma-separated list of thresholds to use", type = str, default = "0.5,0.75" ) 

    args = argp.parse_args()

    thresholds = list( map( float, args.thresholds.split( ',' ) ) )

    sequences = parse_fasta( args.probes )
    moron_scores = scores_to_dict( args.more_ronn )
    iupred_scores = scores_to_dict( args.iupred )

    with open( args.output, 'w' ) as of:
        of.write( 'probe\tmore_ronn_mean\t' + '\t'.join( [ 'more_ronn_cnt_ge_%.2f' % item for item in thresholds ] ) + '\t' +
                  'iupred_mean\t' + '\t'.join( [ 'iupred_cnt_ge_%.2f' % item for item in thresholds ] ) + '\n'
                )

        # for each probe
        for probe in sequences:
            to_write = list()
            mor_list = list()
            iupred_list = list()

            to_write.append( probe.get_name() )
            start, end = probe.get_locations()
            # get the start, end locations
            # get the scores in the range [start, end)
            try:
                mor_score = moron_scores[ probe.name_no_loc() ][ start:end ]
                iupred_score = iupred_scores[ probe.name_no_loc() ][ start:end ]
            except KeyError:
                print( probe.get_name() )


            # calculate mean
            mor_mean = np.mean( mor_score )
            iupred_mean = np.mean( iupred_score )

            mor_list.append( mor_mean )
            iupred_list.append( iupred_mean )

            mor_list += count_positions_above( mor_score, *thresholds )
            iupred_list += count_positions_above( iupred_score, *thresholds )

            to_write += mor_list
            to_write += iupred_list

            of.write( '\t'.join( map( str, to_write ) ) + '\n' )


def count_positions_above( data, *positions ):
    count_ge = lambda data, thresh: sum( [ x >= thresh for x in data ] )
    return [ count_ge( data, c ) for c in positions ]


def scores_to_dict( fname ):
    to_float = lambda x: map( float, x.strip().split( ',' ) )

    out = dict()
    with open( fname, 'r' ) as of:
        for lineno, line in enumerate( of ):
            if lineno: # skip header
                seq_name, scores = line.strip().split( '\t' )
                out[ seq_name ] = list( to_float( scores ) )

    return out

def parse_fasta( fname ):
    names, sequences = oligo.read_fasta_lists( fname )

    return [ SequenceWithLocation( n, s ) for n, s in zip( names, sequences ) ]

class Sequence:
    def __init__( self, name = "", seq = "" ):
        self.name = name
        self.seq = seq

    def get_name( self ):
        return self.name

    def get_seq( self ):
        return self.seq

class SequenceWithLocation( Sequence ):
    def __init__( self, name = "", seq = "" ):
        super().__init__( name = name, seq = seq )

        self.start = 0
        self.end = 0

    def get_locations( self ):
        if not self.start:
            self.start = self.find_start()
        if not self.end:
            self.end = self.find_end()
        return ( self.start, self.end )

    def find_start( self ):
        return int( self.name.split( '_' )[ -2 ] )

    def find_end( self ):
        return int( self.name.split( '_' )[ -1 ] )

    def name_no_loc( self ):
        return '_'.join( self.name.split( '_' )[ :-2 ] )
        

if __name__ == '__main__':
    main()
