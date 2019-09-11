#!/usr/bin/env python3
import argparse
import re
import protein_oligo_library as oligo

def main():
    argp = argparse.ArgumentParser( description = "Given a file noting changes made in ids, "
                                                  "change a map and fasta files."
                                  )

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

    for fasta in fastas:
        fixed_n, fixed_s = fix_fasta_names( fasta, id_map )

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

def is_peptide_name( name ):
    pep_name_re = r'(\S|\s)*OXX=(\S)*_\d+_\d+'
    return re.search( pep_name_re, name )

if __name__ == '__main__':
    main()
