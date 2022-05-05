#!/usr/bin/env Rscript

# Randomly select sites along reference genome
# random samples, random chromosome, random position

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--ref_fai", help="reference genome index file in fasta.fai format")
parser$add_argument("--outfile", help="output file: a table where mutations will be created")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")
parser$add_argument("--sample_size", help="the total number of mutations to be created", default=100000, type='integer')
parser$add_argument("--pos_dist", help="the minimum distance between two mutations in each sample", 
                    default=600, type='integer')

args <- parser$parse_args()
ref_file <- args$ref_fai # ref_file = '/data/home/ywang120/myData/PeterData/reference_genome/dmel_r5.44/dmel-all-chromosome-r5.44.fasta.fai'
ped_file <- args$pedfile # ped_file = '/data/home/ywang120/myData/PeterData/main_pipeline/idata/ped/dmel_fam1.ped'
chr_file <- args$chrfile # chr_file='/data/home/ywang120/myData/PeterData/main_pipeline/idata/call_variants/dmel/fam1/chromosome.config'
out_file <- args$outfile

sample_size = args$sample_size
pos_dist = args$pos_dist

#chr_list = c('2L', '2R', '3L', '3R', 'X')
options(scipen=999) # void scentific notation for pos
#############################

getPos <- function(sampleid, chr, current_dt, chr_len){
    pos <- sample(pos_dist:(chr_len-pos_dist), 1)
    while(TRUE){
        # There should be no created mutation in the same sample and the same region.
        nn <- current_dt[ipos %between% c(pos - pos_dist, pos + pos_dist) &  sid==sampleid & ichr==chr,.N]
        # although the pipeline can detect non-biallelic mutations,
        # we avoid the mutations at the same site in the simulation
        oo = current_dt[ichr == chr & ipos == pos,.N]
        if(nn==0 & oo == 0){
            return(pos)
        }
        pos <- sample(1:chr_len, 1)
    }
}


##############################

ped <- fread(ped_file, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}

chr <- fread(chr_file, header=F, sep=',')
chr_list <- chr[V1 %in% c('AUTO', 'X'), V2]



male_offspring <- ped[V5==1, V2]
female_offspring <- ped[V5==2, V2]
sample_list <- c(male_offspring, female_offspring)

fa <- fread(ref_file, header=F, sep='\t')
fa <- fa[V1 %in% chr_list]
fa[,prop:=V2/sum(V2)]
###########################
# empirical distribution for indel length

# distri_file = 'dmel_fam1_alt_distri_indel.txt'

# dd <- fread(distri_file, header=F, sep='\t')




###########################
randsite <- data.table(indel_type=sample(c('insertion','deletion'), size=sample_size, replace = T),
                       indel_len=sample(1:25, size=sample_size, replace=T), # IMPORTANT: uniform distribution, not using the empirical distribution
                       sid = sample(sample_list, size=sample_size, replace = T),
                       snpidx = 1:sample_size,
                       ichr = sample(fa$V1, size=sample_size, prob=fa$prop, replace = T),
                       ipos = rep(0, sample_size))

for(i in 1:sample_size){
    sampleid <- randsite$sid[i]
    chr <- randsite$ichr[i]
    pos <- getPos(sampleid, chr, randsite, fa[V1==chr,V2])
    randsite[snpidx==i, ipos:=pos]
    randsite[snpidx==i, ichr:=chr]
}

fwrite(randsite, out_file, row.names=F, col.names=T, sep='\t')