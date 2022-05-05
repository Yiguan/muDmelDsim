#!/usr/bin/env python3
import allel
import pandas as pd
import numpy as np
import argparse
import sys

#################################
parser = argparse.ArgumentParser()
parser.add_argument('-invcf',   help='input vcf file that need to be phased')
parser.add_argument('-pedfile', help='pedigree file in PED format describe the relationships')
parser.add_argument('-outfile', help='output file name in compressed format, eg .gz')
parser.add_argument('-chrfile', help='comma-separated chromosome file')
args = parser.parse_args()

## first two samples should be parents, the remaining are offspring!!!!!
# vcf_file = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/recomb/dmel/fam1/dmel_fam1.snp.hardFiltered.biallelic.GQ20.vcf.gz'
# ped_file = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/ped/dmel_fam1.ped'
vcf_file = args.invcf
ped_file = args.pedfile
out_file = args.outfile
chr_file = args.chrfile
#################################


ichr = pd.read_csv(chr_file, header=None)
CHR_LIST = list(ichr[1][ichr[0].isin(['AUTO' ,'X'])])
X_chr = list(ichr[1][ichr[0]=='X'])
autosomes = list(ichr[1][ichr[0]=='AUTO'])


# autosomes = ['2L','2R','3L','3R', '4']
# X_chr = 'X'

#################################
#### parse PED file
ped = pd.read_csv(ped_file, header=None, sep='\t')
motherid = set(ped[3])
if len(motherid)!=1:
    sys.exit('Multiple mother samples inferred from PED file')
motherid = list(motherid)[0]

fatherid = set(ped[2])
if len(fatherid)!=1:
    sys.exit('Multiple mother samples inferred from PED file')
fatherid = list(fatherid)[0]
female_offspring = np.array(ped[1][ped[4]==2])

#################################

aa = allel.read_vcf(vcf_file)
## first two samples should be parents, the remaining are offspring!!!!!
## first is mother, second is father, which we will check later
if not ((motherid in aa['samples'][0:2]) & (fatherid in aa['samples'][0:2])):
    sys.exit("Parents should be the first two samples in the vcf file!")

loc = np.array([x in CHR_LIST for x in aa['variants/CHROM']])
gt = allel.GenotypeArray(aa['calldata/GT'])
gt = gt[loc,:,:]
ac = gt.count_alleles()
loc_seg = ac.is_segregating()
gt_pass = gt.subset(loc_seg)[:]

if aa['samples'][0] == fatherid and aa['samples'][1] == motherid:
    gt_pass[:,0,:], gt_pass[:,1,:] = gt_pass[:,1,:], gt_pass[:,0,:].copy() # swap first two columns to make 'Lady First'
    aa['samples'][0],aa['samples'][1] = aa['samples'][1], aa['samples'][0]


## let's create a psudo X for males inherited from father
## in order to phase X from mother
## allel.phase_by_transmission seems unable to phase haploid
x_loc = np.array([x in X_chr for x in aa['variants/CHROM'][loc][loc_seg]])
x_gt = gt_pass[x_loc]
for i in x_gt[:,1,:]:
    if i[0] == -1:
        i[0] = 0
    i[1] = i[0]

f_hap = x_gt[:,1,1]
for i in range(2,x_gt.shape[1]):
    x_gt[:,i,1] = f_hap

gt_pass[x_loc,:,:] = x_gt


#### Phase 

gt_pass_phased = allel.phase_by_transmission(gt_pass, window_size=100)
gt_pass_phased_str = gt_pass_phased.to_gt().decode('utf-8')

### Output

outdf = pd.DataFrame(data=gt_pass_phased_str, columns=aa['samples'])
outdf['POS'] = aa['variants/POS'][loc][loc_seg]
outdf['CHROM'] = aa['variants/CHROM'][loc][loc_seg]
outdf['REF'] = aa['variants/REF'][loc][loc_seg]  
outdf['ALT'] = aa['variants/ALT'][loc][loc_seg][:,0]
col_order = ['CHROM','POS','REF','ALT']
col_order.extend(aa['samples'])
outdf = outdf.reindex(columns=col_order)
outdf.to_csv(out_file, index=False, compression='gzip')

############################################################
#
# Comment out the following code if need to phase X in females
#
###############################################################
## Phase X chromosome for female offspring
# female_offspring_index = [i for i,j in enumerate(aa['samples']) if j in female_offspring]
# female_offspring_index.insert(0,0)
# female_offspring_index.insert(1,1)

# x_loc = aa['variants/CHROM'] == X_chr
# x_gt = gt[x_loc,:,:]
# x_ac = x_gt.count_alleles()
# x_loc_seg = x_ac.is_segregating()
# x_gt_pass = x_gt.subset(x_loc_seg)[:]
# x_gt_pass_female_offspring = x_gt_pass[:,female_offspring_index,:]
# x_gt_pass_phased = allel.phase_by_transmission(x_gt_pass_female_offspring, window_size=100)

# x_gt_pass_phased_str = x_gt_pass_phased.to_gt().decode('utf-8')

# s_name = [j for i,j in enumerate(aa['samples']) if i in female_offspring_index]
# x_outdf = pd.DataFrame(data=x_gt_pass_phased_str,  columns=s_name)
# x_outdf['POS'] = aa['variants/POS'][x_loc][x_loc_seg]
# x_outdf['CHROM'] = aa['variants/CHROM'][x_loc][x_loc_seg]
# x_outdf['REF'] = aa['variants/REF'][x_loc][x_loc_seg]
# x_outdf['ALT'] = aa['variants/ALT'][x_loc][x_loc_seg][:,0]
# x_col_order = ['CHROM','POS','REF','ALT']
# x_col_order.extend(s_name)
# x_outdf = x_outdf.reindex(columns=x_col_order)
# x_outdf.to_csv("gatk_allsample_scikit_allel_phased.csv.gz", index=False, compression='gzip')
