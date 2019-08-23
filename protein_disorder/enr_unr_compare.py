#!/usr/bin/env python3
import re
import argparse
import scipy.stats
import pandas as pd

def main():
    argp = argparse.ArgumentParser( description = "Compare the " )

    argp.add_argument( '--deconv','-d', help = "Output file produced by PepSIRF's species deconvolution "
                       "module.",
                       type = str
                     )
    argp.add_argument( '--enriched_peptides', '-e', help = "File containing names of significantly enriched peptides, "
                       "one per line.", type = str
                     )
    argp.add_argument( '--peptide_scores', help = "File output by 'summarize_per_probes.py' that contains scores "
                       "for the desired quantiles for each probe.", type = str
                     )
    argp.add_argument( '--pep_name_map', help = "Name of file mapping encoded peptide names to their original name. "
                       "This argument is only necessary if the peptides included in '--enriched_peptides' are encoded. "
                       "Tab-delimited file with the first entry the un-encoded name, and the second the encoded name."
                     )

    args = argp.parse_args()

    peptide_scores = get_peptide_scores( args.peptide_scores )
    enr_species    = pd.read_csv( args.deconv, index_col = 'Species ID', sep = '\t' )
    enr_probes     = set( get_probes_from_file( args.enriched_peptides,
                                                map_fname = args.pep_name_map
                                              )
                        )


    enr_df = pd.DataFrame()


    # get all of the probes, aggregate the counts to a species-level
    species_means_global = aggregate_means( peptide_scores )

    map( lambda x: enr_df.append( species_means_global[ x ] ), enr_probes )
    print( enr_df )

    # get only the enriched probes, aggregate the counts to a species level
    # species_means_enriched = aggregate_means( enr_species )

def get_probes_from_file( probe_fname, map_fname = "" ):
    out = list()
    name_map = dict()


    if map_fname:
        with open( map_fname, 'r' ) as of:
            for lineno, line in enumerate( of ):
                if lineno:
                    spl = line.strip().split( '\t' )
                    name_map[ spl[ 1 ] ] = spl[ 0 ]
                
                       
    with open( probe_fname, 'r' ) as of:
        for line in of:
            if name_map:
                out.append( name_map[ line.strip() ] )
            else:
                out.append( line.strip() )
    return out
            

def get_peptide_scores( fname ):
    return pd.read_csv( fname, sep = '\t' )

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
