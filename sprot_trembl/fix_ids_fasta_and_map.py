#!/usr/bin/env python3
import argparse
import re
import protein_oligo_library as oligo

def main():
    argp = argparse.ArgumentParser( description = "Given a file noting changes made in ids, "
                                                  "change a map and fasta files."
                                  )

    argp.add_argument( '--substitution_record' )
    argp.add_argument( '--fasta', help = "Comma-separated list of fasta files in which to change ids." )
    argp.add_argument( '--map', help = "Name of map file associating "
                                       "long peptide names to short peptide names to change."
                     )
    argp.add_argument( '--suffix', help = "Suffix to add to the end of names of original "
                       "files, before the filetype separator. For example, file.fasta "
                       "becomes file_suffix.fasta",
                       default = "changed_ids"
                     )

    args = argp.parse_args()

    fastas = args.fasta.split( ',' )
    sub_record = parse_sub_record( args.substitution_record )

    for fasta in fastas:
        fixed_n, fixed_s = fix_fasta_names( fasta, id_map )

        prefix, suffix = fasta.split( '.' )
        prefix += '_' + args.suffix 

        oligo.write_fasta_lists( prefix + '.' + suffix, fixed_n, fixed_s )

def parse_sub_record( fname ):
    with open( fname, 'r' ) as of:
        out = dict()
        for lineno, line in enumerate( of ):
            if lineno:
                spl = line.strip().split( '\t' )

                if is_peptide_name( spl[ 0 ] ):
                    # remove location annotation

                    spl[ 0 ] = '_'.join( spl[ 0 ].split( '_' )[ :-2 ] )
                    spl[ 1 ] = '_'.join( spl[ 1 ].split( '_' )[ :-2 ] )
                out[ spl[ 0 ] ] = spl[ 1 ]
        return out

def fix_fasta_names( fname, id_map ):
    name, seqs = oligo.read_fasta_lists( fname )
    fixed_names = list()

    for name in names:
        new_name = name
        if is_peptide_name( name ):
            split = name.split( '_' )
            name_only = '_'.join( split[ :-2 ] )
            location  = '_'.join( split[ -2: ] )

            if name_only in id_map:
                new_name = id_map[ name_only ]
                new_name += location

        else:
            if name in id_map:
                new_name = id_map[ name ]

        fixed_names.append( new_name )
    return fixed_names, seqs

def is_peptide_name( name ):
    pep_name_re = r'(\s|\S)+(_(\d+)_(\d+)){1}'
    return re.search( pep_name_re, name )

if __name__ == '__main__':
    main()
