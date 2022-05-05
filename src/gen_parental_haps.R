#!/usr/bin/env Rscript
library("data.table")
library("dplyr")
library("tidyr")
library("Biostrings")
library("stringr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="input file with phased vcf-liked body, seperated by comma and with header")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--reffa", help="reference genome file in fasta format")
parser$add_argument("--outfile", help="output parental haps in fasta format")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")

args <- parser$parse_args()
infile <- args$infile # infile = 'dmel_fam1.snp.hardFiltered.biallelic.GQ20.phased.gz'
pedfile <- args$pedfile
fafile <- args$reffa
outfile <- args$outfile
chrfile <- args$chrfile


################################################

chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]

chr_list <- c(AUTO_CHR, X_CHR)

species = stringr::str_split(infile, '_')[[1]][1]
fam = paste0('fam', stringr::str_split(infile, '_')[[1]][2])

###################################
# R	        A or G
# Y	        C or T
# S	        G or C
# W	        A or T
# K	        G or T
# M	        A or C
unphased2Het <- function(b1,b2){
    b1 <- toupper(b1)
    b2 <- toupper(b2)
    b12 <- c(b1,b2)
    if(all(b12 %in% c("A","G"))) {ret <- "R"}
    else if(all(b12 %in% c("C","T"))) {ret <- "Y"}
    else if(all(b12 %in% c("G","C"))) {ret <- "S"}
    else if(all(b12 %in% c("A","T"))) {ret <- "W"}
    else if(all(b12 %in% c("G","T"))) {ret <- "K"}
    else if(all(b12 %in% c("A","C"))) {ret <- "M"}
    else{ret <- "N"}
    return(ret)
}

# fasta sequence to data frame
fa2df <- function(chr){
    pp <- paste0("^", chr,"$")
    #seqid <- which(str_detect(stringr::str_split(names(fa),' '), pp))
    tt <- lapply(stringr::str_split(names(fa),' '), function(x) str_detect(x[[1]][1], pp))
    seqid <- which(unlist(tt))
    if(length(seqid)!=1){stop('chromosome name error!')}
    seq <- paste(fa[seqid])
    vdf <- data.table("POS"=1:nchar(seq), "seq" = strsplit(seq,"")[[1]])
    return(vdf)
}

#############################################
ped <- fread(pedfile, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}

aa <- fread(infile, header=T, sep=",")
names(aa) <- stringr::str_replace_all(names(aa), fs, 'father') 
names(aa) <- stringr::str_replace_all(names(aa), ms, 'mother') 



# to fix phasing error in parents
if(!file.exists('tofix.txt')){
    cat('Cannot find file: tofix.txt\n')
    cat('Continue assuming no phasing errors!')
}else{
    fix <- fread('tofix.txt', header = TRUE, sep='\t')
    if(unique(fix$SPECIES)!=species){stop('species name in tofix.txt not match!')}
    if(unique(fix$FAM)!=fam){stop('family name in tofix.txt not match!')}
    for(i in 1:nrow(fix)){
        ichr <- fix[i, CHROM]
        ip <- fix[i, PARENT]
        istart <- fix[i, fix_s]
        iend <- fix[i, fix_e]
        cat('To fix:\n')
        print(fix[i,])
        if(ip == ms){
            aa[CHROM==ichr & POS %between% c(istart, iend), h1:=substr(mother,1,1)]
            aa[CHROM==ichr & POS %between% c(istart, iend), hh:=substr(mother,2,2)]
            aa[CHROM==ichr & POS %between% c(istart, iend), h2:=substr(mother,3,3)]
            aa[CHROM==ichr & POS %between% c(istart, iend), mother:=paste0(h2,hh,h1)]
            aa$h1 <- NULL
            aa$hh <- NULL
            aa$h2 <- NULL
        }else if (ip == fs) {
            aa[CHROM==ichr & POS %between% c(istart, iend), h1:=substr(father,1,1)]
            aa[CHROM==ichr & POS %between% c(istart, iend), hh:=substr(father,2,2)]
            aa[CHROM==ichr & POS %between% c(istart, iend), h2:=substr(father,3,3)]
            aa[CHROM==ichr & POS %between% c(istart, iend), father:=paste0(h2,hh,h1)]
            aa$h1 <- NULL
            aa$hh <- NULL
            aa$h2 <- NULL
        } else{
            stop('patneral ID error!')
        }    
    }
}
names(aa) <- stringr::str_replace_all(names(aa), 'father', fs) 
names(aa) <- stringr::str_replace_all(names(aa), 'mother', ms) 

fwrite(aa, stringr::str_replace(infile, '.gz','_fix.gz'), row.names=F, col.names=T, sep=',', quote=F)


fa <- readDNAStringSet(fafile)


for(chr in chr_list){
    seqdf <- fa2df(chr)
    print(chr)
    for(fm in c(fs, ms)){
        print(fm)
        tmp <- aa %>% select("CHROM","POS","REF","ALT",grep(paste0(fm), names(aa), value = TRUE)) %>%
            filter(CHROM==chr)
        names(tmp)[5] <- c("p.GT")
        # Assume biallelic SNP!!!!!
        tmp <- tmp %>% extract(p.GT,c("h1","phase","h2"),"(.)(.)(.)") 
        # missing(".") to REF, "0 1" to GT
        tmp[h1==".", h1:=REF]
        tmp[h1=="0", h1:=REF]
        tmp[h1=="1", h1:=ALT]
        tmp[h2==".", h2:=REF]
        tmp[h2=="0", h2:=REF]
        tmp[h2=="1", h2:=ALT]
        # unphased het to IPAUC code
        tmp[h1!=h2 & phase=="/", h1:= apply(tmp[h1!=h2 & phase=="/"],1,function(x) unphased2Het(x[3],x[4]))]
        tmp[h1!=h2 & phase=="/", h2:=h1]
        
        tmp_join <- left_join(seqdf, tmp, by = "POS")
        tmp_join[is.na(h1),h1:=seq]
        tmp_join[is.na(h2),h2:=seq]
        # for X in father, it's haploid
        if(chr %in% X_CHR & fm == fs){
            seqinr::write.fasta(sequences = as.list(paste0(tmp_join$h1,collapse="")),
                    names=paste0(chr,"_",fm,"_h1"),
                    nbchar = 50, as.string = TRUE, open="a", 
                    file.out=outfile)
        } else{ 
            seqinr::write.fasta(sequences = as.list(c(paste0(tmp_join$h1,collapse=""), paste0(tmp_join$h2, collapse=""))),
                                names=c(paste0(chr,"_",fm,"_h1"), paste0(chr,"_",fm,"_h2")),
                                nbchar = 50, as.string = TRUE, open="a", 
                                file.out=outfile)
        }
    }
}