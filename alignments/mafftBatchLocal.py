#!/usr/bin/env python

import optparse
from subprocess import Popen, PIPE

# In batch, submits slurm jobs to generate mafft alignments for multiple fasta files

def main():
    usage = '%prog [options] fasta1 [fasta2 ...]'
    p = optparse.OptionParser(usage)
    
    p.add_option('-a', '--alg', default="ginsi,einsi", help='Comma-separated list of algorithms to run. [ginsi,einsi]')
    
    opts, args = p.parse_args()
    
    for i,each in enumerate(args):
        for a in opts.alg.split(","):
            outstr = "%s_%03d" % (a,i)
            cmd="%s --quiet %s >%s_mafft-%s.fasta" % (a, each, ".".join(each.split(".")[:-1]), a)
            print (cmd)
            job=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)

#----------------------End of main()
    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

