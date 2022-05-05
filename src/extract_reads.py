#!/usr/bin/env python3

# select all the reads that contain the mutations in a sample
# output sam file format

import subprocess
import os
import argparse
import sys
import multiprocessing

parser = argparse.ArgumentParser()
parser.add_argument('-infile',  help='a file describing the locations of a site')
parser.add_argument('-bamPath', help='a directory of bam files')

args = parser.parse_args()

in_file = args.infile # in_file = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/denovo_mutations/dmel/fam1/denovo_sim_snp/dmel_fam1_random_site_mutations_snp.txt'
bam_path = args.bamPath # bam_path='/data/home/ywang120/myData/PeterData/main_pipeline/idata/call_variants/dmel/fam1'


# def extractReads(sampleid, chrom, pos):
#     sid = sampleid.split('_')[1]
#     refind = [x for x in os.listdir(bam_path) if x.endswith("_"+sid+".sort.bam")]
#     if len(refind) != 1:
#         sys.exit("Could find bam file or non-unique bam files for sample: {}".format(sampleid))
#     outsam = sampleid+'_'+chrom+'_'+pos+'.sam'
#     fout = open(outsam, 'w')
#     bam = os.path.join(bam_path, refind[0])
#     subprocess.call(['samtools', 'view', bam, chrom+':'+pos+'-'+pos], stdout=fout)
#     fout.close()

def extractReads_multi(itask):
    sampleid = itask[0]
    chrom = itask[1]
    pos = itask[2]
    refind = [x for x in os.listdir(bam_path) if x.endswith("_"+sampleid+".sort.bam")]
    if len(refind) != 1:
        sys.exit("Could find bam file or non-unique bam files for sample: {}".format(sampleid))
    outsam = sampleid+'_'+chrom+'_'+pos+'.sam'
    fout = open(outsam, 'w')
    bam = os.path.join(bam_path, refind[0])
    subprocess.call(['samtools', 'view', bam, chrom+':'+pos+'-'+pos], stdout=fout)
    fout.close()


def main():
    with open(in_file, 'r') as ff:
        lines = ff.readlines()[1:]
        task_list = []
        for ll in lines:
            print(ll)
            lv = ll.strip().split('\t')
            sampleid = lv[0]
            chrome = lv[2]
            pos = lv[3]
            tmp = [sampleid, chrome, pos]
            task_list.append(tmp)
            #extractReads(sampleid, chrome, pos)
    mypool = multiprocessing.Pool(12) # do no use too many cores, I/O too busy
    mypool.map(extractReads_multi,task_list)

if __name__ == "__main__":
    main()
