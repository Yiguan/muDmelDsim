#!/usr/bin/env Rscript


## Summary the distribution of indels given the parents being 
## high confidence homozygous for ref and alt.

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="eg:dmel_fam1.indel.table.gz")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")
parser$add_argument("--outfile", help="output table with alterative read depth")

args <- parser$parse_args()
pedfile <- args$pedfile 
infile <- args$infile
chrfile <- args$chrfile
outfile <- args$outfile
MIN_PARENTAL_DEPTH = 15
MAX_INDEL_LEN = 25 # arbitrary value, the maximum length of INDELs
MAX_DEPTH = 150 # to be summerised 

##################################
ped <- fread(pedfile, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}

male_offspring <- ped[V5==1, V2]
female_offspring <- ped[V5==2, V2]

chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]

##################################
aa <- fread(infile, header = TRUE, sep = "\t")
# remove non-biallelic sites
aa <- aa[!ALT %like% ',']
a_auto <- aa[CHROM %in% AUTO_CHR]

names(a_auto) <- stringr::str_replace_all(names(a_auto), ms, "mother")
names(a_auto) <- stringr::str_replace_all(names(a_auto), fs, "father")

a_auto <- a_auto[mother.NR >= MIN_PARENTAL_DEPTH & father.NR >= MIN_PARENTAL_DEPTH]
a1 <- a_auto[mother.NV==0 & father.NV==father.NR]
a2 <- a_auto[father.NV==0 & mother.NV==mother.NR]
a12 <- rbind(a1,a2)
a12 <- a12[nchar(REF) <= MAX_INDEL_LEN & nchar(ALT)<=MAX_INDEL_LEN]


# as not too many sites, we do not do stratified summary!!
nr <- a12[,paste0(ped$V2,'.NR'), with=F] %>% as.matrix
nv <- a12[,paste0(ped$V2,'.NV'), with=F] %>% as.matrix

outdt <- data.table(dd = 1:MAX_DEPTH, alt='xxx')

for(i in 1:MAX_DEPTH){
    print(i)
    altv <- nv[nr==i]
    print(altv)
    alt_out <- paste(altv, collapse=',')
    outdt[i, alt:=alt_out]
}

fwrite(outdt, outfile, col.names=F, row.names=F, sep='\t', quote=F)