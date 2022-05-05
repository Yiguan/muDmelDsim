#!/usr/bin/env Rscript

library("data.table")
library("dplyr")
library("argparse")
parser <- ArgumentParser()


parser$add_argument("--ref_fai", help="reference genome index file")
parser$add_argument("--mutation_table", help="The original table of all random generated sites")
parser$add_argument("--recover_file", help='Result file after running denovo filtering on simulated data.')
parser$add_argument("--callable_outfile", help="ouput file for callable sites")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="comma separated file, telling whether auto or X chromosome")

args <- parser$parse_args()

fai_file <- args$ref_fai #'/data/home/ywang120/myData/bbsrc/idata/input_dat/ref_genome/dmelanogaster/dmel-all-chromosome-r6.42.fasta.fai'
sim_sites <- args$mutation_table #'dmel_27_random_site_mutations_snp.txt'
sim_result <- args$recover_file #'sim_snp_denovo.txt'
callable_outfile <- args$callable_outfile
chrfile = args$chrfile
pedfile = args$pedfile

###

START_EXTEND = 5 


###
chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]

fai <- fread(fai_file, header=F, sep='\t')
fai <- fai[V1 %in% c(AUTO_CHR, X_CHR)]
##
ped <- fread(pedfile, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}
SAMPLE_N <- nrow(ped)+2
male_offspring <- ped[V5==1, V2]
female_offspring <- ped[V5==2, V2]


#######
fai <- fread(fai_file, header=F, sep='\t', select=1:2)
names(fai) <- c('CHROM', 'LEN')

aa <- fread(sim_sites, header=T, sep='\t')
names(aa) <- c('TYPE', 'INDEL_LEN', 'SID', 'snpidx','CHROM', 'POS')
bb <- fread(sim_result, header=T, sep='\t') 

rname <- bb %>% select(ends_with('NV')) %>% names %>% stringr::str_remove(.,'.NV')

bb$SID <- bb %>% select(ends_with('NV')) %>% apply(.,1, which.max) %>%
    rname[.]
bb <- bb %>% mutate(TYPE = ifelse(nchar(REF)>nchar(ALT),'deletion','insertion')) %>% 
    mutate(LEN=nchar(REF)-nchar(ALT)) %>%
    select(c('TYPE','LEN', 'SID', 'CHROM', 'POS'))  

aa$pos_low <- aa$POS - START_EXTEND
aa$pos_high <- aa$POS + START_EXTEND

aa$called <- 0
bb$simulated <- 0

for(i in 1:nrow(aa)){
    ichr = aa[i,CHROM]
    isid = aa[i, SID]
    ipos_low = aa[i, pos_low]
    ipos_high = aa[i, pos_high]
    if(bb[CHROM==ichr & SID==isid & POS %between% c(ipos_low, ipos_high),.N] >= 1){
        aa[i,called:=1]
        bb[CHROM==ichr & SID==isid & POS %between% c(ipos_low, ipos_high),simulated:=1]
    }
}

sim_nN <- aa[,.(calln=.SD[called==1,.N], simN=.SD[,.N]),by=.(SID, CHROM)] %>%
    mutate(call_prop = calln/simN) %>% 
    left_join(.,fai, by='CHROM') %>%
    mutate(call_len = LEN*call_prop) %>% as.data.table

sim_nN[CHROM %in% AUTO_CHR, total_call_len:=2*call_len]
sim_nN[CHROM %in% X_CHR & SID %in% male_offspring, total_call_len:=call_len]
sim_nN[CHROM %in% X_CHR & SID %in% female_offspring, total_call_len:=2*call_len]

cat('\nThe total callable proportion:\n')
cat(sum(sim_nN$calln)/sum(sim_nN$simN), '\n')

cat('\nThe false discovery rate:\n')
cat(bb[simulated==0,.N]/bb[,.N], '\n')



fwrite(sim_nN, callable_outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')
fwrite(aa, stringr::str_replace(callable_outfile, 'callable', 'checksim'), 
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')

fwrite(bb, stringr::str_replace(callable_outfile, 'callable', 'checkcall'), 
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')