# Scripts used to update the species ids of sequences download from UniProt Sprot/Trembl databases 


### First, locate the changes. Optionally only include changes in sequences found in a fasta.
```
usage: find_changes.py [-h] [--first FIRST] [--second SECOND]
                       [--filter FILTER] [--fasta FASTA] [--output OUTPUT]

Find changes in Taxonomic ID from one version of sprot/tremble data files to
another. Optionally only report results whose ids are in the sequences in a
fasta

optional arguments:
  -h, --help       show this help message and exit
  --first FIRST    The first data file to check.
  --second SECOND  The second data file to check.
  --filter FILTER  Fasta file of sequences to include.
  --fasta FASTA    Optional Fasta file to check for sequences whose ids were
                   changed. If not included, output will contain all of the
                   changed IDs. Otherwise, only the sequences in this file
                   whose ids have changed will be reported.If included, a
                   column will be added to the output with the names of any
                   sequences that were affected.
  --output OUTPUT  Name of file to output.
```

### Update the old ids in a map with new ones.
```
usage: fix_old_ids.py [-h] [--fixed FIXED] [--original ORIGINAL]
                      [--output OUTPUT]

Update old ids in a name map with the new ones.

optional arguments:
  -h, --help           show this help message and exit
  --fixed FIXED        Name of the file that corrects ids (output by
                       find_changes.py)
  --original ORIGINAL  Name of the map file to fix. This file should be a tab-
                       delimited file with two columns. The first column
                       should be the original name, the second the name of the
                       encoded peptide.
  --output OUTPUT      Name of file to output. Format of the output is the
                       same as that of the input, but the names will be fixed
                       to include the correct species-level IDs.
```
