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
    argp.add_argument( '--output', '-o', help = "Name of file to write output to.", default = "species_summaries.tsv" )

    args = argp.parse_args()
    thresholds = list( map( float, args.thresholds.split( ',' ) ) )

    means = pd.read_csv( args.input, sep = '\t' )

    aggreg_means = aggregate_means( means )

    with open( args.output, 'w' ) as of:
        of.write( 'species_id\t' + '\t'.join( [ 'more_ronn_perc_probes_mean_ge_%.2f' %  item for item in thresholds ] )  + '\t' +
                  '\t'.join( [ 'iupred_perc_probes_mean_ge_%.2f' %  item for item in thresholds ] ) + '\n'
        )

        for species, data in aggreg_means.items():
            moron_percents = map( str, percent_items_above( data[ 0 ], *thresholds ) )
            iupred_percents = map( str, percent_items_above( data[ 1 ], *thresholds ) )
            of.write( species + '\t' +  '\t'.join( moron_percents ) + '\t' + '\t'.join( iupred_percents ) + '\n' )

def percent_items_above( data, *thresholds ):
    return [ item / len( data ) for item in count_positions_above( data, *thresholds ) ] 

def count_positions_above( data, *positions ):
    count_ge = lambda data, thresh: sum( [ x >= thresh for x in data ] )
    return [ count_ge( data, c ) for c in positions ]

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
