#!/usr/bin/env Rscript
library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="dmel_27.snp.hardFiltered.table.gz")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="comma separated file, telling whether auto or X chromosome")
parser$add_argument("--maxdepth", help="alt frequency table for depth 1 to maxdepth", default=150, type='integer')
parser$add_argument("--outfile", help="ouput file with alt frequency")

args <- parser$parse_args()
infile <- args$infile
pedfile <- args$pedfile
max_depth <- args$maxdepth
outfile <- args$outfile
chrfile <- args$chrfile
#######################################
getAlt <- function(x){
    alt <- stringr::str_split(x, ',')[[1]][2]
    return(alt)
}

getRef <- function(x){
    ref <- stringr::str_split(x, ',')[[1]][1]
    return(ref)
}

########################################

ped <- fread(pedfile, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}

offspring_list <- ped$V2


chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]

parental_min_dp = 30
parental_max_dp = 200
parental_gq = 70
#########################################


#infile = 'dmel_27.snp.hardFiltered.table.gz'
aa <- fread(infile, header=TRUE, sep='\t')

names(aa) <- stringr::str_replace_all(names(aa), ms, "mother")
names(aa) <- stringr::str_replace_all(names(aa), fs, "father")
aa <- aa[CHROM %in% AUTO_CHR]
# biallelic
aa <- aa[nchar(ALT)==1]
# DP and GQ
aa <- aa[mother.DP %between% c(parental_min_dp, parental_max_dp) &
  father.DP %between% c(parental_min_dp, parental_max_dp) &
  mother.GQ >= parental_gq &
  father.GQ >= parental_gq]
# AD
aa <- aa[(mother.AD %like% '.*,0$' & father.AD %like% '^0,.*') |
   (father.AD %like% '.*,0$' & mother.AD %like% '^0,.*')]

aa<- aa %>% select(paste0(offspring_list,".AD"))

alt_N <- apply(aa, 1:2, getAlt) %>% as.matrix
alt_N <- apply(alt_N, 1:2, as.numeric)

ref_N <- apply(aa, 1:2, getRef) %>% as.matrix
ref_N <- apply(ref_N, 1:2, as.numeric)

dp_N <- alt_N + ref_N

for(dp in 1:max_depth){
    freq_list <- alt_N[dp_N==dp]
    freq_str <- paste(freq_list, collapse=',')
    write(paste(dp, freq_str), outfile, append=TRUE)  
}