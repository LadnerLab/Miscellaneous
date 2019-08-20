#!/usr/bin/env python3
import tempfile
import protein_oligo_library as oligo
import sys
import subprocess


def main():
    if len( sys.argv ) != 2:
        print( "Perform IUPred analysis on a fasta containing "
               "more than one sequence. IUPred must be executable with './iupred2a'. "
               "Output is sent to standard out."
             )

        print( "USAGE: do_pred.py input_file" )

        sys.exit( 1 )

    names, sequences = oligo.read_fasta_lists( sys.argv[ 1 ] )

    for name, sequence in zip( names, sequences ):
        with tempfile.NamedTemporaryFile( mode = 'w' ) as of:
            of.write( '>%s\n%s' % ( name, sequence ) )
            of.flush()

            print( '>' + name )

            print( subprocess.check_output( './iupred2a.py %s long' % of.name, shell = True ).decode( 'ascii' ) )

if __name__ == '__main__':
    main()
