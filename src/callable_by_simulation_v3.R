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
names(aa) <- c('sid', 'snpidx','CHROM', 'POS')
bb <- fread(sim_result, header=T, sep='\t', select=1:5) 


ab <- left_join(aa, bb, by=c('sid','CHROM','POS'))
total_callable <- ab[!is.na(REF),.N]/aa[,.N]
cat('\nThe total callable proportion:\n')
cat(total_callable, '\n')

sim_nN <- ab[!is.na(REF),.(calln=.N),by=.(sid,CHROM)] %>%
    left_join(., ab[,.(simN=.N),by=.(sid, CHROM)], by=c('sid','CHROM')) %>%
    mutate(call_prop = calln/simN) %>% left_join(.,fai, by='CHROM') %>%
    mutate(call_len = LEN*call_prop) %>% as.data.table
names(sim_nN)[1] <-  'SID'
sim_nN[CHROM %in% AUTO_CHR, total_call_len:=2*call_len]
sim_nN[CHROM %in% X_CHR & SID %in% male_offspring, total_call_len:=call_len]
sim_nN[CHROM %in% X_CHR & SID %in% female_offspring, total_call_len:=2*call_len]


aa <- left_join(aa, bb, by=c('sid','CHROM', 'POS')) %>% mutate(called = ifelse(is.na(REF), 0,1)) %>%
    select(-c('REF', 'ALT'))

bb <- left_join(bb, aa, by = c('sid','CHROM', 'POS')) %>% mutate(simulated = ifelse(is.na(snpidx), 0,1)) %>%
    select(-called)

fwrite(sim_nN, callable_outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')

fwrite(aa, stringr::str_replace(callable_outfile, 'callable', 'checksim'), 
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')

fwrite(bb, stringr::str_replace(callable_outfile, 'callable', 'checkcall'), 
       col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')