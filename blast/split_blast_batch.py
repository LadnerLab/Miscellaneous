#! /usr/bin/env python

# By Jason Ladner

# Used to run blasts with multiple queries and/or references
# All combinations of queries and references will be run

from __future__ import division
import sys, optparse, os, math, random
from subprocess import Popen, PIPE
import threading, queue, time, os, subprocess
import inout as io

#Defines a class that will contain all of the info needed to do an individual blsat and parse
class info2blast:
    """A class to store info needed for blast and parse"""
    def __init__(self, opts, this_out, file, cmd):
        prefix='.'.join(file.split('.')[:-1])
        self.reg_out = '%s_parsed.txt' % (prefix)
        self.no_hits = '%s_nohits.txt' % (prefix)
        self.blast_cmd = cmd
        self.parse_cmd = '%s --reg_out %s --no_hits %s --numHits %d --numHsps %d --goodHit %s --xml %s' % (opts.subBlastParse, self.reg_out, self.no_hits, opts.numHits, opts.numHsps, opts.goodHit, this_out)

def main():

    #To parse command line
    usage = "usage: %prog [options] query1 [query2 ...]"
    p = optparse.OptionParser(usage)
    
    #Input/output files
    p.add_option('-s', '--summaryFile', help='Name for summary output file [None, Required]')
    p.add_option('-l', '--logFile', help='Name for log file of commands run [Default will be generated based on --summaryFile]')
    p.add_option('--ns', '--nucSubject', help='Fasta file of nucleotide sequences to compare the query sequences to. Will format if necessary. [None]')
    p.add_option('--ps', '--protSubject', help='Fasta file of protein sequences to compare the query sequences to. Will format if necessary. [None]')
    
    #General
    p.add_option('-n', '--numProcs', type='int', default=4, help='Number of separate blasts to start [4]')
    p.add_option('--temp', default='./temp', help='Name for temporary working directory. This will be created at the begining of the script and then removed at the end [./temp]')
    p.add_option('-b', '--blastType', help='Type of blast to run. Options "blastn", "blastx", "blastp", "tblastx", "tblastn" [blastn or blastx or blastn,blastx]')
    p.add_option('--dontIndex', default=False, action='store_true',  help="Use this flag if you don't want the script to try to index the database. This is necessary for complex databases like nt and nr")
    p.add_option('--keepOut', default=False, action='store_true',  help="Use this flag if you don't want to delete the non-parsed blast files [automatically used if out format not xml]")
    #p.add_option('--dontParse', default=False, action='store_true',  help="Use this flag if you don't want to parse the blast output [automatically used if out format not xml]\n Without parsing it is not possible to do iterative blasting on only non-matching sequences.\n Any iterative blasting will be done on the full query.")
    p.add_option('--blastFull', default=False, action='store_true',  help="Blast full query for each task [automatically used if out format not xml]")
    p.add_option('--blastDBv4', default=False, action='store_true',  help="Includ -blastdb_version 4 when formatting databases [False]")
    p.add_option('--subBlastParse', default="sub_blast_parse_v1.2.py",  help="How to call sub_blast_parse.py ['sub_blast_parse_v1.2.py']")
    
    #Blast options
    p.add_option('--task', default='megablast,dc-megablast,blastn', help='Type of blastn to run. Options "blastn", "dc-megablast", "megablast" [megablast,dc-megablast,blastn]')
    p.add_option('--evalue', default='10', help='Maximum evalue for hit to be recorded [10]')
    p.add_option('--outFmt', type='int', default=5, help="Integer specifying blast output type. 5=xml, 6=tabular, 0=default(I think). Parsing won't work unless xml output is used [5]")
    p.add_option('--numHits', type='int', default=5, help='Integer specifying the number of blast hits to report per query. [5]')
    p.add_option('--numHsps', type='int', default=1, help='Integer specifying the number of alignments to report per query/subject pair. [1]')

    #To determine what goes to next blast stage
    p.add_option('--goodHit', default='0.05', help='Floating point number specifying the evalue nec. for a hit to be significant. [0.05]')
    
    opts, args = p.parse_args()

    #Adjust defaults if output type is not xml
    if opts.outFmt!=5:
        opts.keepOut=True
        opts.dontParse=True
    
    # Set name for log file
    if not opts.logFile:
        opts.logFile = "%s.log" % (opts.summaryFile)
    
    # Open log file for wriitng
    log = open(opts.logFile, "w")
    
    #Read in references
    if opts.ps:
        refD = io.fileDictHeader(opts.ps, "Ref", "Name")
        refType = "protein"
    elif opts.ns:
        refD = io.fileDictHeader(opts.ns, "Ref", "Name")
        refType = "nucleotide"
    
    # Open summary file for writing
    with open(opts.summaryFile, "w") as fout:
        fout.write("Blast Result\tProtein\tQuery Fasta\n")
    
        #Step through each provided reference
        for ref, pname in refD.items():
        
            if refType == "protein":
                opts.ps = os.path.abspath(ref)
            elif refType == "nucleotide":
                opts.ns = os.path.abspath(ref)
        
            #Step through each query
            for query in args:
                opts.query = query
            
                #Change all filenames to their absolute path versions
                opts.query=os.path.abspath(opts.query)
    
                #Save current working directory
                opts.startDir=os.getcwd()

                #Set default blast type if none provided
                if not opts.blastType:
                    if opts.ns and opts.ps: opts.blastType='blastn,blastx'
                    elif opts.ns: opts.blastType='blastn'
                    elif opts.ps: opts.blastType='blastx'
                    else: '!!!ERROR: must provide at least one subject fasta!!!!'
        
                #Step through each type of blast
                for b_type in opts.blastType.split(','):
                    if b_type=='blastn':
                        for bn_type in opts.task.split(','):
                            outName = split_blast(b_type, bn_type, log, opts)
                            #Write out information for this query:reference pair
                            fout.write("%s\t%s\t%s\n" % (outName, pname, query))

                    else:
                        outName = split_blast(b_type, '', log, opts)

                        #Write out information for this query:reference pair
                        fout.write("%s\t%s\t%s\n" % (outName, pname, query))
        
    log.close()

###------------------------End of main()--------------------------------

###--------------------Combine multiple input queries------------------------

def combine_queries(queries, new_name='combo_query_%d.fasta' % random.randrange(9999)):
    fout=open(new_name, 'w')
    for query in queries.split(','):
        fin=open(query, 'r')
        for line in fin:
            fout.write(line)
        fin.close()
    fout.close()
    return new_name


###---------------Split blast function----------------------------------------------------
    
def split_blast (blast_type, task, log, opts):
#    print (blast_type, task)
    #Create temporary working direcory and move to that directory
    if not os.path.exists(opts.temp): os.mkdir(opts.temp)
    os.chdir(opts.temp)
    
    #Create multiple subset query fastas in the temp working directory
    sub_files=split_fasta(opts)
    if sub_files:
        #Check to see is the subject's fasta is formatted as a blast database. If not, format it
        #This section also determines whether the nucleotide or protein reference should be used
        if blast_type in ['blastn', 'tblastx', 'tblastn']:
            subject = opts.ns
            #Check to make sure --dontIndex flag was not used
            if not opts.dontIndex:
                if not os.path.isfile('%s.nsq' % opts.ns):
                    cmd="makeblastdb -in %s -dbtype nucl" % (opts.ns)
                    if opts.blastDBv4:
                        cmd+=" -blastdb_version 4"
                    log.write("%s\n" % cmd)
                    format_db=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                    format_db.wait()
        elif blast_type in ['blastx', 'blastp']:
            subject = opts.ps
            #Check to make sure --dontIndex flag was not used
            if not opts.dontIndex:
                if not os.path.isfile('%s.psq' % opts.ps):    
                    cmd="makeblastdb -in %s -dbtype prot" % (opts.ps)
                    if opts.blastDBv4:
                        cmd+=" -blastdb_version 4"
                    log.write("%s\n" % cmd)
                    format_db=Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
                    format_db.wait()
        #Run and parse blasts
        reg_files=[]
        nohit_files=[]
        blast_result_files=[]
        for file in sub_files: 
            #Only blastn needs the 'task' variable
            if blast_type=='blastn':
                this_out='%s_%s_%s_%s' % (file, blast_type, task[:2], subject.split('/')[-1])
                cmd = '%s -query %s -db %s -evalue %s -out %s -outfmt %d -task %s' % (blast_type, file, subject, opts.evalue, this_out, opts.outFmt, task)
                log.write("%s\n" % cmd)
                blast_result_files.append(this_out)
            else:
                this_out='%s_%s_%s' % (file, blast_type, subject.split('/')[-1])
                cmd = '%s -query %s -db %s -evalue %s -out %s -outfmt %d' % (blast_type, file, subject, opts.evalue, this_out, opts.outFmt)
                log.write("%s\n" % cmd)
                blast_result_files.append(this_out)
                
            work_info = info2blast(opts, this_out, file, cmd)
            request_work(work_info)
            reg_files.append(work_info.reg_out)
            nohit_files.append(work_info.no_hits)
            
        make_and_start_thread_pool(opts.numProcs)
        stop_and_free_thread_pool()
        
        outName, no_good_hits = combine_outputs(blast_type, task, subject, reg_files, nohit_files, opts)
        outName = os.path.abspath(outName)
        
        if not opts.blastFull and not blast_type=="blastp":
            #Make new query for the next round of blasting. Change opts.query to refer to this new file.
            opts.query=subset_fasta(no_good_hits, blast_type, task, opts)
        
        #Leave working directory and clean it
        if not opts.keepOut:
            #Deleting raw blast output and subset query files
            for file in blast_result_files + sub_files + nohit_files: 
                if os.path.isfile(file):
                    os.remove(file)

    os.chdir(opts.startDir)
    if not opts.keepOut: os.rmdir(opts.temp)
    
    return outName


###------------------------Split fasta into multiple files--------------------------------

def split_fasta(opts):
    #Will hold the names of all the files created
    sub_files=[]
    names, seqs = read_fasta_lists(opts.query)
    num_seqs=len(names)
    if num_seqs>=opts.numProcs: sub_size=int(math.ceil(num_seqs/opts.numProcs))                                          #Rounds up so that the first few subsets might have slightly more than the last
    elif num_seqs>0: 
        opts.numProcs=num_seqs
        sub_size=1
    else: return sub_files
    
    for start in range(0, num_seqs, sub_size):
        sub_names=names[start:start+sub_size]
        sub_seqs=seqs[start:start+sub_size]
        new_filename='%d_%d.fasta' % (start+1, start+sub_size)
        sub_files.append(new_filename)
        write_fasta(sub_names, sub_seqs, new_filename)
    return sub_files

###------------------------Parse blast--------------------------

def combine_outputs(blast_type, task, subject, reg_files, nohit_files, opts):
    #Will be for combining subset files
    
    #Generate name for output
    if task:
        outName = '../%s_%s_%s_%s_parsed.txt' % (os.path.basename(opts.query), blast_type, task[:2], subject.split('/')[-1])
    else:
        outName = '../%s_%s_%s_parsed.txt' % (os.path.basename(opts.query), blast_type, subject.split('/')[-1])
    with open(outName, 'w') as out_parse:
        out_parse.write("Query Name\tQuery Length\tSubject Name\tSubject Length\tAlignment Length\tQuery Start\tQuery End\tSubject Start\tSubject End\tHsp Score\tHsp Expect\tHsp Identities\tPercent Match\tNumber_of_gaps\n")
        for f in reg_files:
            fin=open(f, 'r')
            for line in fin:
                out_parse.write(line)
            fin.close()
            os.remove(f)
    
    return outName, read_files_list(nohit_files)


################THREADING#####################################

# Globals (start with a capital letter)
Qin  = queue.Queue()
Qout = queue.Queue()
Qerr = queue.Queue()
Pool = []

def report_error():
    ''' we "report" errors by adding error information to Qerr '''
    Qerr.put(sys.exc_info()[:2])

def get_all_from_queue(Q):
    ''' generator to yield one after the others all items currently
        in the queue Q, without any waiting
    '''
    try:
        while True:
            yield Q.get_nowait()
    except queue.Empty:
        raise StopIteration

def do_work_from_queue():
    ''' the get-some-work, do-some-work main loop of worker threads '''
    while True:
        command, item = Qin.get()       # implicitly stops and waits
        if command == 'stop':
            break
        try:
            # simulated work functionality of a worker thread
            if command == 'process':
                # 'item' is whatever was passed into Qin 
                #This is where I need to add in what I want each thread to do
                
                #Run blast
#                print(item.blast_cmd)
                result1 = run_subproc(item.blast_cmd)
                #Parse blast results
#                print(item.parse_cmd)
                result2 = run_subproc(item.parse_cmd)
                
            else:
#                raise ValueError, 'Unknown command %r' % (command)
                print( 'Unknown command %r' % (command))
        except:
            # unconditional except is right, since we report _all_ errors
            report_error()
        else:
            Qout.put(result1)
            Qout.put(result2)

def make_and_start_thread_pool(number_of_threads_in_pool, daemons=True):
    ''' make a pool of N worker threads, daemonize, and start all of them '''
    for i in range(number_of_threads_in_pool):
         new_thread = threading.Thread(target=do_work_from_queue)
         new_thread.setDaemon(daemons)
         Pool.append(new_thread)
         new_thread.start()

def request_work(data, command='process'):
    ''' work requests are posted as (command, data) pairs to Qin '''
    Qin.put((command, data))

def get_result():
    return Qout.get()     # implicitly stops and waits

def show_all_results():
    for result in get_all_from_queue(Qout):
        print ('Result:', result)

def show_all_errors():
    for etyp, err in get_all_from_queue(Qerr):
        print ('Error:', etyp, err)

def stop_and_free_thread_pool():
    # order is important: first, request all threads to stop...:
    for i in range(len(Pool)):
        request_work(None, 'stop')
    # ...then, wait for each of them to terminate:
    for existing_thread in Pool:
        existing_thread.join()
    # clean up the pool from now-unused thread objects
    del Pool[:]


def run_subproc(cmd):
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print ("Child was terminated by signal", retcode)
        elif retcode >0:
            print ("Child returned", retcode)
#    except OSError, e:
    except:
        print ("Execution failed:", e)
    



###------------------------Read and write files-------------------------

def read_files_list(files):
    l=[]
    for f in files:
        fin=open(f, 'r')
        for line in fin:
            l.append(line.strip())
    return l

# Extracts data from a fasta sequence file. Returns two lists, the first holds the names of the seqs (excluding the '>' symbol), and the second holds the sequences
def read_fasta_lists(file):
    fin = open(file, 'r')
    count=0
    
    names=[]
    seqs=[]
    seq=''
    for line in fin:
        line=line.strip()
        if line and line[0] == '>':                #indicates the name of the sequence
            count+=1
            names.append(line[1:])
            if count>1:
                seqs.append(seq)
            seq=''
        else: seq +=line
    seqs.append(seq)
    
    return names, seqs
    
    

#writes a new fasta file
def write_fasta(names, seqs, new_filename):
    fout=open(new_filename, 'w')
    for i in range(len(names)):
        fout.write(">%s\n%s\n" % (names[i], seqs[i]))
    fout.close()


def recursive_join(list, delimiter="\t"):
    ready_to_join = []
    for index, value in enumerate(list):
        if type(value) == type(()) or type(value) == type([]):
            ready_to_join.append(recursive_join(value))
        elif type(value) == type(1) or type(value) == type(1.0): 
            ready_to_join.append(str(value)) 
        else: 
            ready_to_join.append(value)
            
    joined=delimiter.join(ready_to_join)
    return joined

def subset_fasta(no_good_hits, blast_type, task, opts):    
    names, seqs=read_fasta_lists(opts.query)
    simplenames=[name.split('/t')[0] for name in names]
    subnames=[]
    subseqs=[]
    for i in range(len(seqs)):
        if simplenames[i] in no_good_hits:
            subnames.append(names[i])
            subseqs.append(seqs[i])
#    print (blast_type, task)
    new_query_name = '%s/%s_%s_no_good_hits.fasta' % (opts.startDir, blast_type, task)
    write_fasta(subnames, subseqs, new_query_name)
    return new_query_name    

###---------------------------->>>

if __name__ == "__main__":
    main()  
