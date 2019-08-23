#!/usr/bin/env python3

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
    argp.add_argument( '--original_scores', help = "File output by 'summarize_per_species.py' that contains scores "
                       "for the desired quantiles.", type = str
                     )

    args = argp.parse_args()

if __name__ == '__main__':
    main()
