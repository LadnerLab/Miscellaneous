#!/usr/bin/env python3
import argparse


def main():
    arg_parser = argparse.ArgumentParser( description = '' )

    arg_parser.add_argument( '--input', '-i', help = "Name of the file to accept as input. File should be two tab-delimited columns with "
                             "the second column being comma-delimited float values, each being a positional disorder prediction score for "
                             "the corresponding amino acid.",
                             type = str
    )
    arg_parser.add_argument( '--output', '-o', help = "File to write output to. Output will be two tab-delimited columns with the first as "
                             "sequence names. The second will be a comma-delimited set of ranges (For example, '0-5') for the sequence that "
                             "is above or equal to the given threshold.",
                             default = 'collapsed.tsv'
    )
    arg_parser.add_argument( '--threshold', '-t', help = "The threshold for inclusion. Any ranges above the threshold will be output.",
                             type = float, default = 0.5
                           )

    transform = lambda x: ( x.split( '\t' )[ 0 ],
                            collapse( x.split( '\t' )[ 1 ].split( ',' ), args.threshold )
                          )

    args = arg_parser.parse_args() 

    parsed_data = parse_input( args.input, transform_fn = transform )
    write_output( args.output, parsed_data )


def write_output( out_fname, data ):

    range_to_str = lambda x: '-'.join( [ str( item ) for item in x ] )

    with open( out_fname, 'w' ) as of:
        of.write( 'Sequence Name\tRanges\n' )
        for item in data:
            name = item[ 0 ]
            ranges = item[ 1 ]

            ranges = [ range_to_str( item ) for item in ranges ]

            of.write( '%s\t%s\n' % ( name, ','.join( ranges ) ) )

def collapse( values, threshold ):
    values = [ float( item ) for item in values ]

    ranges_above = list()
    current_range = [ -1, -1 ]

    is_default = lambda x: x == [ -1 , -1 ]


    for idx, item in enumerate( values ):
        if item >= threshold:
            # check if this is the start of a range
            if is_default( current_range ):
                current_range[ 0 ] = idx
            # continuation of a range
            else:
                current_range[ 1 ] = idx
                
        else:
            # end of a range - or continuation of a streak below threshold
            # make sure ranges of length 1 are not included
            if not( is_default( current_range ) ) and current_range[ 0 ] < current_range[ 1 ]:
                ranges_above.append( tuple( current_range.copy() ) )
            current_range = [ -1, -1 ]
    return ranges_above
        

        

def parse_input( in_fname, transform_fn = lambda x: x ):
    output = list()
    with open( in_fname, 'r' ) as in_f:
        for lineno, line in enumerate( in_f ):
            if lineno:
                line = line.strip()
                output.append( transform_fn( line ) )
    return output
                


if __name__ == '__main__':
    main()
