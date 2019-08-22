#!/usr/bin/env python3
import argparse
import pandas as pd
import re

def main():
    argp = argparse.ArgumentParser( description = "Aggregate scores at species level, find %% probes meet threshold."
                                  )

    argp.add_argument( '--input', '-p', help = "File containing probe summaries"
                     )

    argp.add_argument( '--thresholds', '-t', help = "Comma-separated list of thresholds to use", type = str, default = "0.5,0.75" ) 
    argp.add_argument( '--output', '-o', help = "Name of file to write output to.", default = "probe_summaries.tsv" )

    args = argp.parse_args()
    thresholds = list( map( float, args.thresholds.split( ',' ) ) )

    means = pd.read_csv( args.input, sep = '\t' )

    aggreg_means = aggregate_means( means )


def aggregate_means( dataframe ):
    data = dict()

    for index, row in dataframe.iterrows():
        spec_id = get_species_id( row[ 'probe' ] )

        if spec_id not in data:
            data[ spec_id ] = ( list(), list() )
        data[ spec_id ][ 0 ].append( row[ 'more_ronn_mean' ] )
        data[ spec_id ][ 1 ].append( row[ 'iupred_mean' ] )

    return data

def get_species_id( name ):
    ox_re = re.compile( 'OXX=(\d*),(\d*),(\d*),(\d*)' )

    groups = ox_re.search( name ) 

    if groups:
        return groups.group( 2 )
    return None
    


if __name__ == '__main__':
    main()
