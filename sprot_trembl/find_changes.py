#!/usr/bin/env python3
import argparse
import protein_oligo_library as oligo
import re

def main():
    argp = argparse.ArgumentParser( description = 'Find changes in Taxonomic ID from one version of '
                                    'sprot/tremble data files to another. Optionally '
                                    'only report results whose ids are in the sequences in '
                                    'a fasta'
                                  )

    argp.add_argument( '--first', help = "The first data file to check." )
    argp.add_argument( '--second', help = "The second data file to check." )
    argp.add_argument( '--filter', help = "Fasta file of sequences to include. " )
    argp.add_argument( '--fasta', help = "Optional Fasta file to check for sequences whose ids were changed. "
                       "If not included, output will contain all of the changed IDs. \nOtherwise, only the sequences in "
                       "this file whose ids have changed will be reported."
                       "If included, a column will be added to the output with the names of any sequences that were "
                       "affected."
                     )
    argp.add_argument( '--output', help = "Name of file to output.", default = "output.tsv" )

    args = argp.parse_args()

    # parse sequences from the first and second (only need ID and OX)
    first_recs  = parse_records( args.first )
    second_recs = parse_records( args.second )

    # slightly trim search space
    to_search = first_recs.keys() & second_recs.keys()

    # parse sequences from fasta, grab IDs for each
    fasta_records = get_records_from_fasta( args.filter )

    # find those that have changed (ID remains, OX different)
    changed_ids = filter( lambda x: first_recs[ x ] != second_recs[ x ],
                          to_search
                        )
    # filter out IDs that are not from the fasta (if included)
    if fasta_records:
        target_ids = set( changed_ids ) & fasta_records.keys() 
    else:
        target_ids = set( changed_ids )

    affected_peptides = None
    if args.fasta:
        peptides = parse_peptides( args.fasta )
        affected_peptides = find_affected_peptides( peptides, target_ids )

    # write the output
    header = f'id\t{args.first}\t{args.second}'

    if args.fasta:
        header += '\taffected_peptides'

    write_output( args.output,
                  first_recs,
                  second_recs,
                  target_ids,
                  peptides = affected_peptides,
                  header = header
                )

def parse_peptides( fname ):
    names, sequences = oligo.read_fasta_lists( fname )
    return names

def find_affected_peptides( pep_list, ids ):
    id_re = re.compile( "ID=([\S]+)\s[\S]+\sOXX=\d*,(\d*),\d*,\d*" )
    output = dict()

    for pep in pep_list:
        search = re.search( id_re, pep )

        if search:
            id = search.group( 1 )
            if id not in output:
                output[ id ] = list()
            output[ id ].append( pep )
    return output

def write_output( fname, first_records, second_records,
                  target_ids, header = "", peptides = None
                ):

    with open( fname, 'w' ) as of:
        if header:
            of.write( f'{header}\n' )

        for id in target_ids:
            of.write( f'{id}\t{first_records[ id ]}\t{second_records[ id ]}' )

            if peptides and id in peptides:
                # join on '~' because that character is not in any names
                of.write( f"\t{ '~'.join( peptides[ id ] ) }" )
            of.write( '\n' )
            

def get_records_from_fasta( fname ):
    records = dict()
    line_re = re.compile( "ID=([\S]+)\s[\S]+\sOXX=\d*,(\d*),\d*,\d*" )
    if fname:
        names, seqs = oligo.read_fasta_lists( fname )
        for name in names:
            name_search = re.search( line_re, name )
            if name_search:
                id = name_search.group( 1 )
                ox = name_search.group( 2 )

                records[ id ] = str( ox )
            else:
                print( "No match: ", name )
    return records

def parse_records( fname ):
    records = dict()

    new_record = lambda x: x[ 0 ] == 'ID'
    record_id  = lambda x: x[ 0 ] == 'OX'
    get_id     = lambda x: x.split( '=' )[ 1 ].replace( ';', '' )

    with open( fname, 'r' ) as of:
        for lineno, line in enumerate( of ):
            spl = line.strip().split()
            try:
                
                if new_record( spl ):
                    records[ spl[ 1 ] ] = ''
                    last_id = spl[ 1 ]
                elif record_id( spl ):
                    records[ last_id ] = get_id( spl[ 1 ] )
            except IndexError:
                print( "No id found for line %d: " % ( lineno + 1 ),  line.strip() )
                
    return records

if __name__ == '__main__':
    main()
