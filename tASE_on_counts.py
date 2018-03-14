#! /usr/bin/env python


'''
tASE_on_counts.py [OPTIONA] g_listfile c_listfile

tASE takes two files with a list of positional counts (from nuc_by_pos.py) for genomic DNA and cDNA
assumes paired g/cDNA. Files should be a list of files, relative paths allowed in combination with -w (see below)


        all 4 below arguments must be used, and multiple samples can be combined, using a , to separate
    -p [INT] position - position to determine ASE over (1-indexed)
    -m [nuc] marine nucleotide
    -f [nuc] fresh nucleotide
    -l [string] locus name

    -o [fpath]  outfile, saved in tsv format
	-d [fpath] data file, storing SL, TL, GT data for each animal
	-w [directory] working directory of the count files
'''

import sys,getopt,module
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Arial')
matplotlib.rc('text', usetex='false')
matplotlib.rc('pdf', fonttype=42)

args = sys.argv[1:]
optlist, args =  getopt.getopt(args, 'p:f:m:l:o:w:d:')


POSITION = 0
LOCUS = ''
FRESH = 'A'
MARINE = 'T'
G_FILE = args[0]
C_FILE = args[1]
OUTFILE = ''
WORKING_DIRECTORY = ''
DATAFILE = ''

NTS = ['A','C','G','T']
NT_INDEX = {'A':0,'C':1,'G':2,'T':3}

for a in optlist:
    if a[0] == '-p':
        POSITION = int(a[1])
    if a[0] == '-f':
        FRESH = a[1]
    if a[0] == '-m':
        MARINE = a[1]
    if a[0] == '-l':
        LOCUS = a[1]
    if a[0] == '-o':
        OUTFILE = a[1]
    if a[0] == '-w':
        WORKING_DIRECTORY = a[1]
    if a[0] == '-d':
        DATAFILE = a[1]


##loads counts from files in file_list as position and nucs (F,M) default
def load_counts(file_list,position,F,M):
    
    ##stores counts of reads in each rep F,M list
    count_list = [[],[]]
    
    flh = open(file_list,'r')
    
    for count_file in flh:
        f = open(WORKING_DIRECTORY + count_file.strip(),'r')
        ##print f
        lines = f.readlines()
        
        ##print lines
        
        F_counts_split = lines[NT_INDEX[F]].split('\t')
        M_counts_split = lines[NT_INDEX[M]].split('\t')
        
        F_counts = float(F_counts_split[position])
        M_counts = float(M_counts_split[position])
        
        count_list[0].append(F_counts)
        count_list[1].append(M_counts)
        f.close()

    return count_list

##plot ratio of the two counts gDNA and cDNA
##show mean and std
def plot_counts(g_counts,c_counts,PLOT=False):
    
    g_xs = np.zeros(len(g_counts[0])) + 1
    c_xs = np.zeros(len(g_counts[0])) + 2
    

    
    g_F_counts = np.array(g_counts[0])
    g_M_counts = np.array(g_counts[1])
    g_ratio = g_F_counts/(g_F_counts + g_M_counts)
    
    
    c_F_counts = np.array(c_counts[0])
    c_M_counts = np.array(c_counts[1])
    c_ratio = c_F_counts/(c_F_counts + c_M_counts)
    
    gmask = np.isfinite(g_ratio)
    cmask = np.isfinite(c_ratio)
    ##print g_ratio,c_ratio,gmask,cmask
    g_ratio = g_ratio[gmask]
    c_ratio = c_ratio[cmask]
    ##print g_ratio,c_ratio
    corrected_g = g_ratio / g_ratio.mean()
    corrected_c = c_ratio / g_ratio.mean()
    
    ##diff= g_ratio / c_ratio
    ##diff = np.log2(g_ratio) - np.log2(c_ratio)
    ##print list(g_ratio),list(c_ratio)
    
    corrected_g = g_ratio / g_ratio.mean()
    corrected_c = c_ratio / g_ratio.mean()
    ##print list(corrected_g),list(corrected_c)
    '''
    print 'g_ratio shaprio:', stats.shapiro(g_ratio)
    print 'g_ratio shaprio,log2:',stats.shapiro(np.log2(g_ratio))
    print 'g_ratio shaprio,arcsin:',stats.shapiro(np.arcsin(g_ratio))
    print 'c_ratio shaprio:',stats.shapiro(c_ratio)
    print 'c_ratio shaprio,log2:',stats.shapiro(np.log2(c_ratio))
    print 'c_ratio shaprio,arcsin:',stats.shapiro(np.arcsin(c_ratio))
    '''
    
    print 't_test:', stats.ttest_ind(g_ratio,c_ratio)
    print 'corrected t_test:', stats.ttest_ind(corrected_g,corrected_c)
    ##print 'paired t_test:', stats.ttest_rel(g_ratio,c_ratio)
    print 't_test,log2:',stats.ttest_ind(np.log2(g_ratio),np.log2(c_ratio))
    print 't_test,differences:',stats.ttest_1samp(corrected_c,0)
    ##print 'Wilcoxon:',stats.wilcoxon(g_ratio,c_ratio)
    
    ##print list(diff)
    
    ymean = [g_ratio.mean(),c_ratio.mean()]
    yerr = [g_ratio.std(),c_ratio.std()]
    if PLOT:
        plt.figure()
        plt.scatter(g_xs,g_ratio,color='blue')
        plt.scatter(c_xs,c_ratio,color='red')
        plt.errorbar([1,2],ymean,yerr=yerr,fmt='.',color='black')
    return g_ratio,c_ratio
    ##plt.show()

##gets sample names from each file and ids
def load_sample_names(infile):
    
    sample_ids = []
    sample_names = []
    
    f = open(infile,'r')
    
    for line in f:
        line = line.strip()
        split = line.split('_')
        cur_name = split[0]
        ##cur_id = split[0][1:]
        cur_id = '_'.join(split[0:3])
        sample_names.append(cur_name)
        sample_ids.append(cur_id)
        
    return sample_names,sample_ids

##laods data from the given file
##stores data in a dictionary, with keys as gDNA and values as (genotype,SL,TL,FISH)
def load_datafile(dfile,locus):
    f = open(dfile,'r')
    
    name_dict = {}
    
    header = f.readline().strip()
    header_list = header.split(',')
    gdna_index = header_list.index(locus + '_gDNA')
    for line in f:
        line = line.strip()
        split = line.split(',')
        gDNA_name = split[gdna_index]
        SL = split[1]
        TL = split[2]
        GT = split[3]
        FISH = split[0]
        name_dict[gDNA_name] = (GT,SL,TL,FISH)
        
    return name_dict

if OUTFILE:
    w = open(OUTFILE,'w')
    w.write('\t'.join(['Sample_Name','Fish','Locus','gDNA','cDNA','corrected_gDNA_ratio','corrected_cDNA_ratio']) + '\n')

if OUTFILE and DATAFILE:
    w = open(OUTFILE,'w')
    w.write('\t'.join(['Sample_Name','Fish','Locus','gDNA','cDNA','corrected_gDNA_ratio','corrected_cDNA_ratio','SL','TL','Genotype']) + '\n')

##print 'san',sample_names,fish_ids
g_counts = load_counts(G_FILE,POSITION,FRESH,MARINE)
c_counts = load_counts(C_FILE,POSITION,FRESH,MARINE)
g_ratio,c_ratio = plot_counts(g_counts,c_counts)
if DATAFILE:
    name_dict = load_datafile(DATAFILE,LOCUS)
    sample_type,sample_ids = load_sample_names(G_FILE)


    #gmask = np.isfinite(g_ratio)
    #cmask = np.isfinite(c_ratio)
    #print g_ratio,c_ratio,gmask,cmask
    #g_ratio = g_ratio[gmask]
    #c_ratio = c_ratio[cmask]
    #print g_ratio,c_ratio

    corrected_g = g_ratio / g_ratio.mean()
    corrected_c = c_ratio / g_ratio.mean()

    for indI,i in enumerate(sample_ids):
        gt,sl,tl,fish_name = name_dict[i]
        write_str = '\t'.join(map(lambda x: str(x),
        [i,fish_name,LOCUS,g_ratio[indI],c_ratio[indI],corrected_g[indI],corrected_c[indI],sl,tl,gt]))
        if OUTFILE:
            w.write(write_str + '\n')
'''

##stores the fc_information for each fish
fish_fc = {}
##for each file mentioned
overall_corrected_c = []
for gfile,cfile,m_base,f_base,position,locus in zip(G_FILE,C_FILE,MARINE,FRESH,POSITION,LOCUS):
    
    sample_names,fish_ids = load_sample_names(gfile)
    print 'san',sample_names,fish_ids
    g_counts = load_counts(gfile,position,f_base,m_base)
    c_counts = load_counts(cfile,position,f_base,m_base)
    g_ratio,c_ratio = plot_counts(g_counts,c_counts)

    #gmask = np.isfinite(g_ratio)
    #cmask = np.isfinite(c_ratio)
    #print g_ratio,c_ratio,gmask,cmask
    #g_ratio = g_ratio[gmask]
    #c_ratio = c_ratio[cmask]
    #print g_ratio,c_ratio

    corrected_g = g_ratio / g_ratio.mean()
    corrected_c = c_ratio / g_ratio.mean()
    overall_corrected_c.append(corrected_c)
    print g_ratio,g_counts
    print sample_names,fish_ids
    print corrected_c
    
    for i in range(len(sample_names)):
        print i
        cur_fish  = fish_ids[i]
        if cur_fish not in fish_fc:
            fish_fc[cur_fish] = []
        fish_fc[cur_fish].append(corrected_c[i])
        write_str = '\t'.join(map(lambda x: str(x),
        [sample_names[i],fish_ids[i],locus,g_ratio[i],c_ratio[i],corrected_g[i],corrected_c[i]]))
        if OUTFILE:
            w.write(write_str + '\n')
    
##G_COUNTS = load_counts(G_FILE,POSITION,FRESH,MARINE)
##C_COUNTS = load_counts(C_FILE,POSITION,FRESH,MARINE)
print overall_corrected_c
print fish_fc

##plot the combined cDNA ratios
plt.figure()
for i,locus in enumerate(LOCUS):
    corrected_c = overall_corrected_c[i]
    xs = np.zeros(len(corrected_c)) + i + 1
    plt.scatter(xs,corrected_c)
    plt.errorbar(i+1,corrected_c.mean(),yerr=corrected_c.std(),fmt='.',color='black')

xs =  np.arange(len(LOCUS)) + 1
plt.xticks(xs,LOCUS)
for fish in fish_fc:
    corrected_c =fish_fc[fish]
    xs =  np.arange(len(corrected_c )) + 1
    plt.plot(xs,corrected_c,color='black')
    
if OUTFILE:
    plot_out = OUTFILE.split('.')[:-1]
    plot_out = '.'.join(plot_out) + '.pdf'
    plt.savefig(plot_out)
plt.show()

'''
    