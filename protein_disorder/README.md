### transform_output
```
usage: transform_output.py [-h] [--input INPUT] [--output OUTPUT]

Transform output friom IUPred and moreRONN into a format that is more easily-
parseable.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        The input file to transform. This should be a file
                        either output from IUPred or moreRONN. It will be
                        inferred from which program the output is.
  --output OUTPUT, -o OUTPUT
                        Name of file to write output to. File will contain a
                        header followed by one line per sequence. Each line
                        will be a comma-separated list of disorder predictions
                        cores, one per amino acid in the sequence.
```

### do_pred
```
Perform IUPred analysis on a fasta containing more than one sequence. 
IUPred must be executable with './iupred2a'. Output is sent to standard out.
USAGE: do_pred.py input_file
```

### collapse_ranges
```
usage: collapse_ranges.py [-h] [--input INPUT] [--output OUTPUT]
                          [--threshold THRESHOLD]

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT, -i INPUT
                        Name of the file to accept as input. File should be
                        two tab-delimited columns with the second column being
                        comma-delimited float values, each being a positional
                        disorder prediction score for the corresponding amino
                        acid.
  --output OUTPUT, -o OUTPUT
                        File to write output to. Output will be two tab-
                        delimited columns with the first as sequence names.
                        The second will be a comma-delimited set of ranges
                        (For example, '0-5') for the sequence that is above or
                        equal to the given threshold.
  --threshold THRESHOLD, -t THRESHOLD
                        The threshold for inclusion. Any ranges above the
                        threshold will be output.
```
