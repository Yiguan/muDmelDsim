#!/usr/bin/env python3

# change read base to create mutations
# number of reads to be changed are based firstly on haplotypes and then the empirical distribution of heterzygotes.
# all males X chromosomes are regarded as homo - all reads to be changed

from itertools import starmap
import collections
import gzip
import pandas as pd
import re
import random
import subprocess
import sys
import os
import multiprocessing
import argparse

pd.set_option('mode.chained_assignment', None)

#################################
parser = argparse.ArgumentParser()
parser.add_argument('-ref',   help='the reference genome in fasta format')
parser.add_argument('-pedfile', help='pedigree file in PED format describe the relationships')
parser.add_argument('-vcfgz',   help='vcf file in gz format as SNP markers to phase')
parser.add_argument('-altdistri', help='ALT count in heterozygotes for different depth, each line: DEPTH C1,C2,C1,C3,C2,....')
parser.add_argument('-logfile', help='tab-separated file that records if the mutation successfully created')
parser.add_argument('-chrfile', help='comma-separated chromosome file')
#parser.add_argument('-outfile', help='')

args = parser.parse_args()

PED_FILE = args.pedfile # ped_file = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/ped/dmel_fam1.ped'
REF_FILE = args.ref # REF_FILE = '/data/home/ywang120/myData/PeterData/reference_genome/dmel_r5.44_edit/dmel-all-chromosome-r5.44_edit.fasta'
INPUT_SNP = args.vcfgz # INPUT_SNP = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/call_variants/dmel/fam1/dmel_fam1.snp.hardFiltered.vcf.gz'
ALT_FREQ = args.altdistri # ALT_FREQ = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/denovo_mutations/dmel/fam1/denovo_sim/dmel_fam1_alt_distri.txt'
LOG_OUT = args.logfile
CHR_FILE = args.chrfile
# the following values should match denovo filtering values
# for male offspring, use half of these values
MAX_DEPTH = 150 
MIN_DEPTH = 10 

# CHR_LIST = ['2L', '2R','3L', '3R','X']
# X_CHR = ['X']

###################################
#### parse PED file
ped = pd.read_csv(PED_FILE, header=None, sep='\t')
motherid = set(ped[3])
if len(motherid)!=1:
    sys.exit('Multiple mother samples inferred from PED file')
motherid = list(motherid)[0]

fatherid = set(ped[2])
if len(fatherid)!=1:
    sys.exit('Multiple mother samples inferred from PED file')
fatherid = list(fatherid)[0]
male_offspring = list(ped[1][ped[4]==1])
female_offspring = list(ped[1][ped[4]==2])

### parse CHR file
ichr = pd.read_csv(CHR_FILE, header=None)
CHR_LIST = list(ichr[1][ichr[0].isin(['AUTO' ,'X'])])
X_CHR = list(ichr[1][ichr[0]=='X'])
#####################################

def gzvcf2pd(filename):
    chrl = []
    posl = []
    refl = []
    altl = []
    polymophic_samples = []
    with gzip.open(filename, 'rt', encoding='utf-8') as ff:
        lines = ff.readlines()
        for ll in lines:
            if ll[0:2] == '##':
                continue
            if ll[0:2] == '#C': #lines[58]
                hh = ll.strip().split('\t') 
                continue
            lv = ll.strip().split('\t') # lines[2200]
            if len(lv[4]) != 1: # only use biallelic snps
                continue
            if lv[0] in CHR_LIST:
                chrl.append(lv[0])
                posl.append(lv[1])
                refl.append(lv[3])
                altl.append(lv[4])
                ff = [ii for ii,fi in enumerate(lv) if fi[0:3]=='0/1' \
                    or fi[0:3]=='1/0' or fi[0:3]=='0|1' or fi[0:3]=='1|0']
                snp_sample = [hh[ii] for ii in ff]
                if lv[0] in X_CHR: # remove male samples if genotyped as het on X
                    snp_sample = [x for x in snp_sample if x not in male_offspring]
                polymophic_samples.append(snp_sample)
    posl = [int(x) for x in posl]
    ret = pd.DataFrame({'CHROM':chrl, 'POS':posl, 'REF':refl, 'ALT':altl, 'POLYSAMPLE':polymophic_samples})
    return ret

def getRefBase(chrom, pos):
    ''' get a reference base from REF_FILE(global) using samtools'''
    pos=str(pos)
    refbase=subprocess.run(["samtools", "faidx", REF_FILE, chrom+":"+ pos+"-"+pos], stdout=subprocess.PIPE,universal_newlines=True)
    ref=re.findall(r"\n(.)\n", refbase.stdout)
    return(ref[0].upper())

def genAlt(ref):
    '''randomly select a base except reference base'''
    bases = ["A", "G", "C", "T"]
    bases.remove(ref)
    return(random.sample(bases,1)[0])


def refCount(depth):
    '''determine a number of psudo mutation for a depth'''
    ref_pool = []
    with open(ALT_FREQ, 'r') as ff:
        lines = ff.readlines()
        for ll in lines:
            rr = int(ll.split(' ')[0])
            if rr == depth:
                tt = ll.strip().split(' ')
                if len(tt)==1: #no empirical distribution for alt
                    ref_pool = [int(depth/2)]
                else:
                    ref_pool = tt[1].split(',')
                break
    if len(ref_pool) == 0:
        sys.exit(f'Depth {depth} could not found in empirical distribution!')
    refcount = random.sample(ref_pool,1)[0]
    return(refcount)


def getReadPos(cigar, target_pos, start_pos):
    '''return the position of base (start from 0) on READ '''
    # remove soft clip, which isn't regarded as read start position in SAM file
    # but present in sequence
    # c_letters = re.findall('[A-Z]', cigar)
    # if 'D' not in c_letters and 'I' not in c_letters:
    #     return(target_pos - start_pos)
    c_item = re.findall('[0-9]*[A-Z]', cigar)
    c_list = [[int(x[:-1]),x[-1]] for x in c_item]
    t_list = []
    last_pos = start_pos - 1
    for ii in c_list:
        if ii[1] == 'S':
            pass
        elif ii[1] == 'H':
            pass
        elif ii[1] == 'I':
            t_list.extend([0]*ii[0])
        elif ii[1] == 'D':
            tmp = [x + last_pos + 1 for x in list(range(ii[0]))]
            last_pos = tmp[-1]
        else:
            tmp = [x + last_pos + 1 for x in list(range(ii[0]))]
            last_pos = tmp[-1]
            t_list.extend(tmp)
    try:
        read_index = t_list.index(target_pos)
    except:
        return(-1)
    # offset soft clip
    if c_list[0][1] == 'S':
        read_index = read_index + c_list[0][0]
    return read_index


def phase2SNPs(target_snp, snp1_pos, snp2_pos, samdf):
    t_sam = samdf[(samdf[3]<=snp1_pos) & (samdf['end_pos']>=snp2_pos)]
    if len(t_sam)==0:
        # phase break: not info available to infer phase
        # we randomly assign the phase to connect the breakpoint
        ref1 = target_snp.loc[target_snp['POS'] == snp1_pos, 'REF'].values[0]
        ref2 = target_snp.loc[target_snp['POS'] == snp2_pos, 'REF'].values[0]
        alt1 = target_snp.loc[target_snp['POS'] == snp1_pos, 'ALT'].values[0]
        alt2 = target_snp.loc[target_snp['POS'] == snp2_pos, 'ALT'].values[0]
        snp1 = [ref1, alt1]
        snp2 = [ref2, alt2]
        t1 = random.sample(snp1, 1)[0]
        t2 = random.sample(snp2, 1)[0]
        hap1 = t1 + t2
        snp1.remove(t1)
        snp2.remove(t2)
        hap2 = snp1[0] + snp2[0]
        ret_dict = {hap1:[[0],[snp1_pos, snp2_pos]], hap2:[[0],[snp1_pos,snp2_pos]]}
        return ret_dict
    t_sam = t_sam.reset_index(drop=True)
    #return(t_sam)
    hap_df = []
    for i in range(0,len(t_sam)):
        seq = t_sam[9][i]
        cigar = t_sam[5][i]
        start_pos = t_sam[3][i]
        base1 = seq[getReadPos(cigar, snp1_pos, start_pos)]
        base2 = seq[getReadPos(cigar, snp2_pos, start_pos)]
        tt = base1+base2
        hap_df.append(tt)
    aa = collections.Counter(hap_df)
    if len(aa) == 1:  # in case one hap can be supported by reads
        phased_base = list([x for x in aa.keys()][0])
        tt = zip(target_snp['REF'], target_snp['ALT'])
        tt = [*tt]
        second_phase = []
        for k in range(0, len(phased_base)):
            mm = [m for m in tt[k] if m != phased_base[k]]
            second_phase.append(mm[0])
        aa[''.join(second_phase)] = 1  # Assume one read support the phase
    ret_dict = {k: [[v], [snp1_pos, snp2_pos]] for k, v in aa.items()}
    return ret_dict


def catPhase(dict1, dict2):
    '''
    dict1 = {'CA': [[10], [9906981, 9907041]], 'TT': [[2], [9906981, 9907041]]}
    dict2 = {'AT': [[11], [9907041, 9907102]], 'TN': [[1], [9907041, 9907102]], 'TC': [[6], [9907041, 9907102]]}
    '''
    ret_phase = {}
    for k1 in dict1.keys():
        anchor_base = k1[-1]
        next_phase = {k2:v2 for k2,v2 in dict2.items() if k2[0]==anchor_base}
        if next_phase == {}:
            ret_phase[k1] = dict1[k1]
        #most_likely_phase = max(next_phase, key=next_phase.get) # only using the most likely phase
        for nn in next_phase.keys():
            phased = k1+nn[1]
            support_reads_count = [*dict1[k1][0], *next_phase[nn][0]]
            snp_pos_list = [*dict1[k1][1], next_phase[nn][1][1]]
            ret_phase[phased] = [support_reads_count, snp_pos_list]            
    return ret_phase


def getReadRange(cigar, start_pos):
    '''get the end of read on reference genome based on starting and cigar'''
    c_item = re.findall('[0-9]*[A-Z]', cigar)
    c_list = [[int(x[:-1]),x[-1]] for x in c_item]
    end_pos = start_pos + sum([l for l,c in c_list if c!='I' and c!='S' and c!='H'])
    return(end_pos-1)

def editCharInString(ss, pos, cc):
    new_s = list(ss)
    new_s[pos] = cc
    return("".join(new_s))


def chooseReads(hap_dict, samdf):
    '''hap_dict = {'CATT': [[10, 11, 9], [9906981, 9907041, 9907102, 9907123]]}
    return a list of index for the reads'''
    selected_read = []
    no_phase_read = []
    hap_snp_base = list(hap_dict.keys())[0]
    hap_snp_pos = hap_dict[hap_snp_base][1]
    for i in samdf.index:
        read_start = samdf[3][i]
        read_end = samdf['end_pos'][i]
        cigar = samdf[5][i]
        seq = samdf[9][i]
        hit = False # check if the read has at least 1 snp
        for j in range(0, len(hap_snp_pos)):
            if hap_snp_pos[j] >= read_start and hap_snp_pos[j] <= read_end:
                hit = True
                snp_base = seq[getReadPos(cigar, hap_snp_pos[j], read_start)]
                if snp_base == hap_snp_base[j]:
                    selected_read.append(i)
                    break # if multiple snp markers on a read, using the first one 
        if not hit:
            no_phase_read.append(i)
        ## For read, phase cannot be determined, select based on empirical distribution
    no_phase_depth = len(no_phase_read)
    if no_phase_depth != 0:
        ref_size = refCount(no_phase_depth)
        tt = random.sample(no_phase_read, int(ref_size))
        selected_read.extend(tt)
    return selected_read


def changeBase(sam_file, read_index, samdf, mutate_pos, mutate_base):
    outname = sam_file.replace('.sam','_mutated.sam')
    for i in read_index:
        seq = samdf[9][i]
        start_pos = samdf[3][i]
        cigar = samdf[5][i]
        readpos = getReadPos(cigar, int(mutate_pos), int(start_pos))
        if readpos == -1:
            return -1
        try:
            seq_edited = editCharInString(seq, readpos, mutate_base)
        except:
            print("Error at:\n")
            print(sam_file) 
            print(read_index)
            sys.exit(0)
        samdf[9][i] = seq_edited
    outsam = samdf.drop('end_pos',axis=1)
    outsam.to_csv(outname, index=False, header=False, sep = '\t')
    return 0


def valideDiploid(hap1_dict, hap2_dict):
    '''
    hap1_dict = {'ACTACG': [[17, 19, 14, 12, 8], [6018229, 6018257, 6018264, 6018286, 6018294, 6018303]]}
    hap2_dict = {'ATACTA': [[14, 16, 11, 10, 9], [6018229, 6018257, 6018264, 6018286, 6018294, 6018303]]}
    '''
    hap1 = list(hap1_dict.keys())[0]
    hap2 = list(hap2_dict.keys())[0]
    h12 = [*zip(hap1,hap2)]
    valid_index = [ii for ii,jj in enumerate(h12) if jj[0]!=jj[1]]
    hap1_new = [list(hap1)[i] for i in valid_index]
    hap2_new = [list(hap2)[i] for i in valid_index]
    hap1_snp_pos_new = [hap1_dict[hap1][1][i] for i in valid_index]
    hap2_snp_pos_new = [hap2_dict[hap2][1][i] for i in valid_index]
    # supporting reads count should not be used in the future, assign 0
    ret_list = [{''.join(hap1_new):[[0],hap1_snp_pos_new]}, {''.join(hap2_new):[[0], hap2_snp_pos_new]}]
    return ret_list



def sub_main(SAM_FILE):
    print(SAM_FILE, file=sys.stderr)
    SID = SAM_FILE.split('_')[0]
    CHR_ID = SAM_FILE.split('_')[1]
    POS = SAM_FILE.split('_')[2].replace(".sam","")
    if os.stat(SAM_FILE).st_size == 0:
        return (SAM_FILE,'Fail')
    SAM_DF = pd.read_csv(SAM_FILE, delimiter='\t', header=None, usecols=range(0,11))
    READ_DEPTH = len(SAM_DF)
    SAM_DF['end_pos'] = list(starmap(getReadRange, zip(SAM_DF[5],SAM_DF[3])))
    REF = getRefBase(CHR_ID, POS)
    #######################
    # for male X chromosomes
    if SID in male_offspring and CHR_ID in X_CHR:
        if READ_DEPTH > MAX_DEPTH/2 or READ_DEPTH < MIN_DEPTH/2:
            return (SAM_FILE,'Fail')
        else:
            selected_read_index = list(range(0, READ_DEPTH)) # all reads in male X-chr to be mutated
            psudo_alt = genAlt(REF)
            ret_code = changeBase(SAM_FILE, selected_read_index, SAM_DF.copy(), POS, psudo_alt)
            if ret_code == 0:
                return (SAM_FILE, 'Succeed')
            else:
                return (SAM_FILE, 'Fail') # in deletions
    else:
        if READ_DEPTH > MAX_DEPTH or READ_DEPTH < MIN_DEPTH:
            return (SAM_FILE,'Fail')
        else:
            # subset SNPs within the region
            min_start = min(SAM_DF[3])
            max_end = max(SAM_DF['end_pos'])
            target_snp = SNP_DF[(SNP_DF['POS'] >= min_start) & (SNP_DF['POS'] <= max_end) & (SNP_DF['CHROM']==CHR_ID)]
            target_snp = target_snp.sort_values('POS')
            target_snp = target_snp.reset_index(drop=True)
            target_snp = target_snp[[SID in bb for bb in target_snp['POLYSAMPLE']]] # make sure the selected sample are polymorphic for the SNP sites
            if len(target_snp) == 0:
                ref_size = refCount(READ_DEPTH)
                selected_read_index = random.sample(range(0, READ_DEPTH), int(ref_size))
                psudo_alt = genAlt(REF)
                ret_code = changeBase(SAM_FILE, selected_read_index, SAM_DF.copy(), POS, psudo_alt)
                if ret_code == 0:
                    return (SAM_FILE, 'Succeed')
                else:
                    return (SAM_FILE, 'Fail') # in deletions
            elif len(target_snp) == 1:
                target_snp = target_snp.sort_values('POS')
                target_snp = target_snp.reset_index(drop=True)
                if int(POS) in target_snp.POS.values:  # synetic mutation hit a polymorphic site, 
                    return (SAM_FILE, 'Fail')          # polymorphic sites are regarded as uncallable sites
                snp_ref = target_snp['REF'][0]
                snp_alt = target_snp['ALT'][0]
                snp_pos = target_snp['POS'][0]
                # random select hap
                hap_base = random.sample([snp_ref, snp_alt],1)[0]
                target_hap_dict = {hap_base: [[0], [snp_pos]]}
                selected_read_index = chooseReads(target_hap_dict, SAM_DF)
                psudo_alt = genAlt(REF)
                ret_code = changeBase(SAM_FILE, selected_read_index, SAM_DF.copy(), POS, psudo_alt)
                if ret_code == 0:
                    return (SAM_FILE, 'Succeed')
                else:
                    return (SAM_FILE, 'Fail') # in deletions
            else:
                target_snp = target_snp.sort_values('POS')
                target_snp = target_snp.reset_index(drop=True)
                if int(POS) in target_snp.POS.values:  # synetic mutation hit a polymorphic site
                    return (SAM_FILE, 'Fail')          # polymorphic sites are regarded as uncallable sites
                marker_count  = len(target_snp)
                if marker_count>10:
                    tt = random.sample(range(marker_count),10) # to avoid 2^N memory problem, random select 10 markers to phase
                    target_snp = target_snp.iloc[tt]
                    target_snp = target_snp.sort_values('POS')
                    target_snp = target_snp.reset_index(drop=True)
                for i in range(0, len(target_snp)-1):
                    hap_dict = phase2SNPs(target_snp, target_snp['POS'][i], target_snp['POS'][i+1], SAM_DF)
                    # first snp set
                    if i==0:
                        main_dict = hap_dict
                    else:
                        main_dict = catPhase(main_dict, hap_dict)
                #print(main_dict)
                # Only use the first two haps with most supporting reads
                first_hap = max(main_dict, key = lambda x:sum(main_dict[x][0]))
                first_hap_dict = {k:v for k,v in main_dict.items() if k == first_hap}
                main_dict.pop(first_hap)
                second_hap = max(main_dict, key = lambda x:sum(main_dict[x][0]))
                second_hap_dict = {k:v for k,v in main_dict.items() if k == second_hap}
                # validate diploid, remove non-polymophic SNP markers for diploid
                valid_hap_list = valideDiploid(first_hap_dict, second_hap_dict) #Caution: if two hap different number of snp markers, may cause bugs!
                # randomly select one of the two haps to create mutations
                target_hap_dict = random.sample(valid_hap_list,1)[0]
                selected_read_index = chooseReads(target_hap_dict, SAM_DF)
                psudo_alt = genAlt(REF)
                ret_code = changeBase(SAM_FILE, selected_read_index, SAM_DF.copy(), POS, psudo_alt)
                if ret_code == 0:
                    return (SAM_FILE, 'Succeed')
                else:
                    return (SAM_FILE, 'Fail') # in deletions


def main():
    sams = [x for x in os.listdir('.') if re.match(r".*_[0-9]*.sam$", x)]
    # sam_done = [x for x in os.listdir() if x.endswith('mutated.sam')]
    # sams = [x for x in sams if x.replace('.sam', '_mutated.sam') not in sam_done]
    global SNP_DF
    SNP_DF = gzvcf2pd(INPUT_SNP)
    # for ss in sams:
    #     print(ss, file=sys.stderr)
    #     ret = sub_main(ss)
    #     ff = [x[0] for x in ret]
    #     rr = [x[1] for x in ret]
    #     out = pd.DataFrame({'file':ff, 'result':rr})
    #     out.to_csv(LOG_OUT, encoding = 'utf-8', header='true', index=False, sep='\t')
    pool = multiprocessing.Pool(48)
    ret = pool.map(sub_main, sams)
    ff = [x[0] for x in ret]
    rr = [x[1] for x in ret]
    out = pd.DataFrame({'file':ff, 'result':rr})
    out.to_csv(LOG_OUT, encoding = 'utf-8', header='true', index=False, sep='\t')

if __name__ == "__main__":
    main()
