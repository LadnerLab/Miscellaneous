#!/usr/bin/env python

import optparse
from subprocess import Popen, PIPE

# In batch, submits slurm jobs to generate mafft alignments for multiple fasta files

def main():
    usage = '%prog [options] fasta1 [fasta2 ...]'
    p = optparse.OptionParser(usage)
    
    p.add_option('-a', '--alg', default="ginsi,einsi", help='Comma-separated list of algorithms to run. [ginsi,einsi]')
    p.add_option('-t', '--time', default="6:00", help='String to use to set slurm time limit on each job. Default is for 6 hours. [6:00]')
    p.add_option('-c', '--threads', default=2, type='int', help='Number of threads to use per job. [2]')
    p.add_option('-m', '--mem', default="1G", help='Amount of memory to request for each job. [1G]')
    
    opts, args = p.parse_args()
    
    for i,each in enumerate(args):
        for a in opts.alg.split(","):
            outstr = "%s_%03d" % (a,i)
            cmd="sbatch -J %s -o %s.out -e %s.err -c %d --mem=%s --wrap='%s --quiet --thread %s %s >%s_mafft-%s.fasta'" % (outstr, outstr, outstr, opts.threads, opts.mem, a, opts.threads, each, ".".join(each.split(".")[:-1]), a)
            print (cmd)
            job=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            jobid_ir = job.stdout.read().split()[3]

#----------------------End of main()
    
###------------------------------------->>>>    

if __name__ == "__main__":
    main()

