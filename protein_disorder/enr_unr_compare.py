#!/usr/bin/env python3
import re
import argparse
import scipy.stats
import pandas as pd
import numpy as np

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
    argp.add_argument( '--quantiles', '-q', help = "The quantiles to compute for each moreRONN and IUPred scores. Should "
                       "be comma-delimited list of values.",
                       type = str, default = "0.25,0.5,0.75"
                     )
    argp.add_argument( '--output', '-o', help = "Name of file to write output to.",
                       type = str, default = "enr_unr_output.tsv"
                     ) 

    args = argp.parse_args()

    quantiles = sorted( map( float, args.quantiles.strip().split( ',' ) ) )

    peptide_scores = get_peptide_scores( args.peptide_scores )
    enr_species    = pd.read_csv( args.deconv, index_col = 'Species ID', sep = '\t' )
    enr_probes     = get_probes_from_file( args.enriched_peptides,
                                                map_fname = args.pep_name_map
                                              )
                     
    # get all of the probes, aggregate the counts to a species-level
    species_means_global = aggregate_means( peptide_scores )

    enr_probe_data = get_enr_probe_data( peptide_scores, enr_probes )
    species_means_enr = aggregate_means( enr_probe_data )

    # remove any unenriched species from our consideration
    # use list() for eager evaluation
    list( map( lambda x: species_means_enr.pop( x ),
               filter( lambda x: int( x ) not in enr_species.index, list( species_means_enr.keys() ) )
             )
        )

    quantiles_enr = get_quantiles( species_means_enr, quantiles )
    quantiles_unr = get_quantiles( species_means_global, quantiles )

    p_values = get_p_vals( species_means_global, species_means_enr )

    write_outputs( args.output,
                   q_enr = quantiles_enr,
                   q_unr = quantiles_unr,
                   means_unr = species_means_global,
                   means_enr = species_means_enr,
                   p_scores = p_values,
                   quantiles = quantiles
                 )

def get_p_vals( unenriched, enriched ):
    output = dict()

    for key, means in enriched.items():
        output[ key ] = [ 0, 0 ]

        # moron 
        output[ key ][ 0 ] = scipy.stats.ttest_ind( means[ 0 ], unenriched[ key ][ 0 ] ).pvalue

        #iupred 
        output[ key ][ 1 ] = scipy.stats.ttest_ind( means[ 1 ], unenriched[ key ][ 1 ] ).pvalue
    return output


def write_outputs( of_name,
                   q_enr, q_unr,
                   means_unr, means_enr,
                   quantiles,
                   p_scores
                 ):
    with open( of_name, 'w' ) as of:
        of.write( 'species_id\t' )

        list( map(
                  lambda x: of.write( 'more_ronn_enriched_quantile_%.2f\tiupred_enriched_quantile_%.2f\t' % ( x, x ) ),
                  quantiles
                 )
            )
        list( map(
                  lambda x: of.write( 'more_ronn_unenriched_quantile_%.2f\tiupred_unenriched_quantile_%.2f\t' % ( x, x ) ),
                  quantiles
                 )
            )

        of.write( 'more_ronn_enriched_mean\tiupred_enriched_mean\t' )
        of.write( 'more_ronn_unenriched_mean\tiupred_unenriched_mean\t' )
        of.write( 'more_ronn_p_value\tiupred_p_value\n' )

        for key, quantile in q_enr.items():
            of.write( key + '\t' )
            # write enriched quantiles
            for q_id, q_val in quantile.items():
                of.write( '%.6f\t%.6f\t' % ( q_val[ 0 ], q_val[ 1 ] ) )

            # write unenriched quantiles
            for q_id, q_val in q_unr[ key ].items():
                of.write( '%.6f\t%.6f\t' % ( q_val[ 0 ], q_val[ 1 ] ) )

            # write enriched means
            of.write( '%.6f\t%.6f\t' % ( np.mean( means_enr[ key ][ 0 ] ),
                                         np.mean( means_enr[ key ][ 1 ] )
                                       )
            )

            # write unenriched means
            of.write( '%.6f\t%.6f\t' % ( np.mean( means_unr[ key ][ 0 ] ),
                                         np.mean( means_unr[ key ][ 1 ] )
                                       )
            )


            # write moron and iupred p values
            of.write( '%.6f\t%.6f' % ( p_scores[ key ][ 0 ], p_scores[ key ][ 1 ] ) )
                

            of.write( '\n' )
            


def get_quantiles( data, quantiles ):
    output = dict()
    for id, scores in data.items():
        output[ id ] = dict()

        for q in quantiles:
            if q not in output[ id ]:
                output[ id ][ q ] = [ 0, 0 ]
            output[ id ][ q ][ 0 ] = np.quantile( scores[ 0 ], q )
            output[ id ][ q ][ 1 ] = np.quantile( scores[ 1 ], q )

    return output

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
   
def get_enr_probe_data( df, to_get ):
    ret_df = pd.DataFrame( data = None,
                           columns = df.columns,
                           index = df.index
                         )

    to_get_set = set( to_get )

    for index, row in df.iterrows():
        if row[ 'probe' ] in to_get_set:
            ret_df = ret_df.append( row )
        
    return ret_df.dropna()

if __name__ == '__main__':
    main()
