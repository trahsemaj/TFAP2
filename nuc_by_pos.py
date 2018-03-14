#! /usr/bin/env python
import sys,getopt,module
import matplotlib.pyplot as plt
import numpy as np
from  Bio import pairwise2
'''
nuc_by_pos.py [options] -o counts.out in.fq

takes a in.fq file, outputs the nt counts per postion to counts
in the form: tsv, with rows = nts, columns = positions

    -q [INT] min qual score (30)
    -s [FILE] fasta file with the reference sequence, must have > m proportion matching nts
    -m [FLOAT] score cutoff for pairwise alignment, % of total length (default .9)
    -o [FILE] outfile [default test.tsv]
    -c [INT] cuts the reads to the given length
    -t [INT] trim off first t bases of a read (read will stay length c)
'''
args = sys.argv[1:]
optlist, args =  getopt.getopt(args, 'q:s:m:o:c:t:')

QUAL_CUTOFF = 30
MATCH_CUTOFF = .95
SEQ_FILE = ''
REF_SEQS = []
OUT_FILE = 'test.tsv'
READ_FILE = args[0]
TRIM = 0



READ_CUT = 46
len_reads = 0
NTS = ['A','C','G','T']

for a in optlist:
    
    if a[0] == '-q':
        QUAL_CUTOFF = int(a[1])
    if a[0] == '-m':
        MATCH_CUTOFF = float(a[1])
    if a[0] == '-s':
        SEQ_FILE = a[1]
    if a[0] == '-o':
        OUT_FILE = a[1]
    if a[0] == '-c':
        READ_CUT = int(a[1])
    if a[0] == '-t':
        TRIM = int(a[1])
        
##load the ref sequence if specified
if SEQ_FILE:
    seq_dict = module.loadRef(SEQ_FILE)

    ref_ids = seq_dict.keys()
    REF_SEQS = []
    for ref in ref_ids:
        REF_SEQS.append(seq_dict[ref][:READ_CUT])

##takes two strings of the same length
##returns the proportion identity between the two (by pairwise matching, no indels)
def pairwise_match(string1,string2):
    matching = 0
    total = len(string2)
    
    for c1,c2 in zip(string1,string2):
        if c1 == c2:
            matching +=1
    
    ##print string1,string2, float(matching)/total
    return float(matching)/total


##takes two strings of the same length
##returns the proportion identity between the two (by pairwise alignment)
def pairwise_align(string1,string2):

    alignments = pairwise2.align.globalxx(string1,string2)
    alignments.sort(key=lambda x:int(x[2]),reverse=True)
    best_alignment = alignments[0]
    best_score = best_alignment[2]
    
    return float(best_score)/len(string1)


##counts the bases within a given read
##outputs two lists - a freq and a count
def count_reads():
    ##dictionary which stores nt_counts
    nt_positions = {}
    total_reads = 0
    filtered_reads = 0
    filtered_bases= 0
    counted_reads = 0
    counted_bases= 0
    ##Read in the read file 4 lines at a time
    f = open(READ_FILE,'r')
    nt_seq = ' '
    while nt_seq:
        REV_COMP = False
        header = f.readline().strip()
        nt_seq = f.readline().strip()[TRIM:READ_CUT+TRIM]
        f.readline()
        qual_scores = f.readline().strip()[TRIM:READ_CUT+TRIM]
        
        ##initialize the dictionary
        if not nt_positions:
            len_reads = len(nt_seq)
            ##print len_reads
            for i in range(len(nt_seq)):
                nt_positions[i] = {}
                for nuc in NTS:
                    nt_positions[i][nuc] = 0
        
        total_reads += 1
        ##print REF_SEQS
        if REF_SEQS:
            fw_cur_matching =  pairwise_match(nt_seq,REF_SEQS[0])
            m_cur_matching = pairwise_match(nt_seq,REF_SEQS[1])
            ##print fw_cur_matching,m_cur_matching
            ##rv_cur_matching = pairwise_match(module.revComp(nt_seq),REF_SEQ)
            if max(fw_cur_matching,m_cur_matching) < MATCH_CUTOFF:
                filtered_reads += 1
                continue
            
            
        
        ##print fw_cur_matching, nt_seq
        
        counted_reads += 1
        
        ##for each nuc
        for cur_index,cur_qual_char,cur_nt in zip(range(len(qual_scores)),qual_scores,nt_seq):
            ##print cur_index,cur_nt, cur_qual_char, ord(cur_qual_char), ord(cur_qual_char)-33
            cur_qual_score = ord(cur_qual_char)-33
            if cur_qual_score < QUAL_CUTOFF:
                filtered_bases += 1
                continue
            
            nt_positions[cur_index][cur_nt] += 1
            counted_bases += 1


    print >> sys.stderr, 'Total Reads: %s' % total_reads
    print >> sys.stderr, 'Filtered Reads: %s' % filtered_reads
    print >> sys.stderr, 'Counted Reads: %s' % counted_reads
    print >> sys.stderr, 'Filtered Bases: %s' % filtered_bases
    print >> sys.stderr, 'Counted Bases: %s' % counted_bases

    ##generate an outlist

    ##A,C,G,T
    
    nt_pos_list = [[],[],[],[]]
    nt_freq_list = [[],[],[],[]]

    for i in range(len_reads):
        cur_nt_counts = nt_positions[i]

        cur_counted = sum(cur_nt_counts.values())
        for index,nt in enumerate(NTS):
            nt_pos_list[index].append(cur_nt_counts[nt])
            if cur_counted > 0:
                nt_freq_list[index].append(float(cur_nt_counts[nt])/cur_counted)
            else:
                nt_freq_list[index].append(0.)
            
    return nt_pos_list,nt_freq_list


##writes the given counts/freq list to a file
def write_counts(OUTFILE,list):

    w = open(OUTFILE,'w')



    for index,count in enumerate(list):
        count = count[:]
        count.insert(0,NTS[index])
        count = map(lambda x:str(x),count)
        towrite = '\t'.join(count)
        w.write(towrite + '\n')
        

    w.flush()
    w.close()



        
def plot_freq(freq_list,OUTFILE=''):
    
        
        
    ##a a row name
    xs = range(len(freq_list[0]))
    for i in range(4):
        
        
        plt.plot(xs,freq_list[i],label=NTS[i])
        
        

    plt.xlim([0,len(freq_list[0])])
    plt.legend()
    plt.grid()
    plt.savefig(OUT_FILE + '.pdf')
    ##plt.show()

nt_counts,nt_freq = count_reads()
write_counts(OUT_FILE+ '_counts.tsv',nt_counts[:])
write_counts(OUT_FILE+ '_freq.tsv',nt_freq[:])
plot_freq(nt_freq,OUTFILE=OUT_FILE + '.pdf')
    




