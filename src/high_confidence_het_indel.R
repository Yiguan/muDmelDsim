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
MIN_PARENTAL_DEPTH = 30
MAX_INDEL_LEN = 25 # arbitrary value, the maximum length of INDELs
MAX_DEPTH = 150

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

###############################

aa <- fread(infile, header = TRUE, sep = "\t")
# remove non-biallelic sites
aa <- aa[!ALT %like% ',']
a_auto <- aa[CHROM %in% AUTO_CHR]

names(a_auto) <- stringr::str_replace_all(names(a_auto), ms, "mother")
names(a_auto) <- stringr::str_replace_all(names(a_auto), fs, "father")

a_auto$mother.GT <- stringr::str_replace_all(a_auto$mother.GT, '\\|','/')
a_auto$father.GT <- stringr::str_replace_all(a_auto$father.GT, '\\|','/')

a_auto_parentalHomo <- a_auto[(mother.GT==paste0(REF,'/',REF) & father.GT==paste0(ALT,'/',ALT)) | 
        (father.GT==paste0(REF,'/',REF) & mother.GT==paste0(ALT,'/',ALT))]
a_auto_parentalHomo <- a_auto_parentalHomo[mother.DP >= MIN_PARENTAL_DEPTH &
                                               father.DP >= MIN_PARENTAL_DEPTH]

a_auto_parentalHomo[,indel_len := nchar(REF) - nchar(ALT)]
a_auto_parentalHomo <- a_auto_parentalHomo[indel_len %between% c(-MAX_INDEL_LEN, MAX_INDEL_LEN)]

ad <- a_auto_parentalHomo[,.SD,.SDcols=c('CHROM', 'POS', 'indel_len',
                                         grep('^s.*.AD$', names(a_auto_parentalHomo), value = TRUE))]
ad <- ad %>% tidyr::gather(sid, val, grep('^s.*.AD$', names(ad), value = TRUE)) %>%
    tidyr::separate( val, into=c('ref','alt'), sep=',', convert = TRUE) %>% as.data.table
ad[,abslen:=abs(indel_len)]
ad[,depth:=ref+alt]
ad$indel_type <- ifelse(ad$indel_len>0,'deletion','insertion')
ad <- ad[depth<=MAX_DEPTH]

ad[abslen<=5, len_bin:='(0-5]']
ad[abslen > 5 & abslen<=10, len_bin:='(5-10]']
ad[abslen > 10 & abslen<=15, len_bin:='(10-15]']
ad[abslen > 15 & abslen<=20, len_bin:='(15-20]']
ad[abslen > 20 & abslen<=25, len_bin:='(20-25]']
ad_id <- unique(ad[,c('indel_type','len_bin','depth')])
for(i in 1:nrow(ad_id)){
    t1 <- ad_id[[i,'indel_type']]
    t2 <- ad_id[[i,'len_bin']]
    t3 <- ad_id[[i, 'depth']]
    tmp <- ad[indel_type==t1 & len_bin==t2 & depth==t3]
    ad_id[i, alt_sum:=paste0(tmp$alt, collapse = ',')]
}

fwrite(ad_id, outfile, col.names=F, row.names=F, sep='\t', quote=F)