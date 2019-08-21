#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import csv
import math
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# IMPORT MY LATEX SO I CAN USE \TEXTSC
import matplotlib as mpl
#mpl.rc('text', **{'usetex':True})
import re

def parse_values( f1, f2 ):
    transf = lambda x: [ float( item ) for item in x.split( ',' ) ]
    x, y = list(), list()
    with open( f1, 'r' ) as of, open( f2, 'r' ) as of2:
        count = 0
        for line1, line2 in zip( of, of2 ):
            if count:
                spl1 = line1.split( '\t' )
                spl2 = line2.split( '\t' )
                x += transf( spl1[ 1 ] )
                y += transf( spl2[ 1 ] )
            count += 1
    return x, y
def main():
    iupred_in = 'iupred_output_transformed.tsv'
    moron_in  = 'moron_output_transformed.tsv'

    x, y = parse_values( iupred_in, moron_in )
    fig = plt.figure( figsize = ( 5, 4 ) )
    ax = fig.add_subplot( 111 )
    ax.scatter( x, y )
    ax.set_xlabel( 'IUPred Scores')
    ax.set_ylabel( 'MoreRONN Scores')

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]
    
    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    plt.savefig( 'correlation.pdf', bbox_inches = 'tight' )

if __name__ == '__main__':
    main()
