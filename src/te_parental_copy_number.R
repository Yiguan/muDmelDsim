#!/usr/bin/env Rscript

# Searching for TEe in parental samples and count copy number in parents.
# homozygous as two copies, heterozygous as one copy

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="Teflon output gt file")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--outfile", help="ouput file with alt frequency")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")
parser$add_argument("--AB", help="TE balance for each loci")
args <- parser$parse_args()
infile <- args$infile # infile='/data/home/ywang120/myData/PeterData/main_pipeline/idata/TE/dmel/fam1/genotypes/dmel_fam1_te_gt.txt'
pedfile <- args$pedfile # pedfile= '/data/home/ywang120/myData/PeterData/main_pipeline/idata/ped/dmel_fam1.ped'
outfile <- args$outfile
chrfile <- args$chrfile
MIN_AB <- as.numeric(args$AB)

MIN_TOTAL_READ = 10
MIN_READ = 2
#MIN_AB = 0.2  ## [0.2-0.8] is considered as het
################################################

chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]
Y_CHR <- chr[V1=='Y', V2]


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

offspring_list <- ped$V2
male_samples <- c(fs, ped[V5==1, V2])
female_samples <- c(ms, ped[V5==2, V2])

sample_N <- nrow(ped) + 2

############################################
aa <- fread(infile, header=F, sep="\t")
a_auto <- aa[V1 %in% AUTO_CHR]

idf <- a_auto[,c(1:4,14)] %>% unique

a_auto_te <- data.table()

for(i in 1:nrow(idf)){
    #print(i)
    chr = idf[i,V1]
    p5 = idf[i, V2]
    p3 = idf[i,V3]
    tid = idf[i,V4]
    tuid = idf[i, V14]
    tmp <- a_auto[V1==chr & V2==p5 &  V3==p3 &V4==tid & V14==tuid]
    # a few duplicate lines
    if(tmp[,.N]!=sample_N){
        if(tmp[,.N]!=sample_N*2){
            warning("strange number of samples! Skipped...")
            next
        }else{
            tmp <- tmp[V7 %like% ","] # duplication caused by 'V7'
        }
    }
    if(tmp[,.N]==0) next
    # ALL the samples have a sufficient reads
    if(!all(tmp[,V10+V11] >= MIN_TOTAL_READ)) next
    if(any(tmp[,V13]==-9)) next
    # the TE must present in either father or mother or both
    if(!(tmp[V15==fs,V10] >= MIN_READ | tmp[V15==ms,V10] >= MIN_READ)) next
    # maternal counts
    if(tmp[V15==ms,V13] %between% c(MIN_AB, 1-MIN_AB)) {
        tmp$te_source <- 'maternal_het'
        a_auto_te <- rbind(a_auto_te, tmp[V15==ms])
    } else if(tmp[V15==ms,V13] > (1-MIN_AB)) {
        tmp$te_source <- 'maternal_hom'
        a_auto_te <- rbind(a_auto_te, tmp[V15==ms])
    }else {}
    # paternal counts
    if(tmp[V15==fs,V13] %between% c(MIN_AB, 1-MIN_AB)) {
        tmp$te_source <- 'paternal_het'
        a_auto_te <- rbind(a_auto_te,  tmp[V15==fs])
    } else if(tmp[V15==fs,V13] > (1-MIN_AB)) {
        tmp$te_source <- 'paternal_hom'
        a_auto_te <- rbind(a_auto_te,  tmp[V15==fs])
    } else{}

}

##################################################
ax <- aa[V1 %in% X_CHR]

idf <- ax[,c(1:4,14)] %>% unique

ax_te <- data.table()
for(i in 1:nrow(idf)){
    #print(i)
    chr = idf[i,V1]
    p5 = idf[i, V2]
    p3 = idf[i,V3]
    tid = idf[i,V4]
    tuid = idf[i, V14]
    tmp <- ax[V1==chr & V2==p5 &  V3==p3 &V4==tid & V14==tuid]
    # a few duplicate lines
    if(tmp[,.N]!=sample_N){
        if(tmp[,.N]!=sample_N*2){
            warning("strange number of samples! Skipped...")
            next
        }else{
            tmp <- tmp[V7 %like% ","] # duplication caused by 'V7'
        }
    }
    if(tmp[,.N]==0) next
    # all females should meet minimum read depth on X
    if(!all(tmp[V15 %in% female_samples,V10+V11] >= MIN_TOTAL_READ)) next
    # all males should meet minimum read depth/2 on X
    if(!all(tmp[V15 %in% male_samples,V10+V11] >= MIN_TOTAL_READ/2)) next
    if(any(tmp[,V13]==-9)) next
    # the TE must present in either father or mother or both
    if(!(tmp[V15==fs,V10] >= MIN_READ | tmp[V15==ms,V10] >= MIN_READ)) next
    # maternal counts
    if(tmp[V15==ms,V13] %between% c(MIN_AB, 1-MIN_AB)) {
        tmp$te_source <- 'maternal_het'
        ax_te <- rbind(ax_te, tmp[V15==ms])
    } else if(tmp[V15==ms,V13] > (1-MIN_AB)) {
        tmp$te_source <- 'maternal_hom'
        ax_te <- rbind(ax_te, tmp[V15==ms])
    }else {}
    # paternal counts
    if(tmp[V15==fs,V13] >= (1-MIN_AB)) {
        tmp$te_source <- 'paternal_het' # using this to indict one copy on X
        ax_te <- rbind(ax_te, tmp[V15==fs])
    }  else{}
}

######################################################

ay <- aa[V1 %in% Y_CHR]

ay_te <- data.table()
if(ay[,.N]!=0){
    idf <- ay[,c(1:4,14)] %>% unique
    for(i in 1:nrow(idf)){
        #print(i)
        chr = idf[i,V1]
        p5 = idf[i, V2]
        p3 = idf[i,V3]
        tid = idf[i,V4]
        tuid = idf[i, V14]
        tmp <- ay[V1==chr & V2==p5 &  V3==p3 &V4==tid & V14==tuid]
        # a few duplicate lines
        if(tmp[,.N]!=sample_N){
            if(tmp[,.N]!=sample_N*2){
                warning("strange number of samples! Skipped...")
                next
            }else{
                tmp <- tmp[V7 %like% ","] # duplication caused by 'V7'
            }
        }
        if(tmp[,.N]==0) next
        # all males should meet minimum read depth/2 on X
        if(!all(tmp[V15 %in% male_samples,V10+V11] >= MIN_TOTAL_READ/2)) next
        if(any(tmp[V15%in% male_samples, V13]==-9)) next
        # the TE must present in  father 
        if(!(tmp[V15==fs,V10] >= MIN_READ)) next
        # paternal counts
        if(tmp[V15==fs,V13] >= (1-MIN_AB)) {
            tmp$te_source <- 'paternal_het' # using this to indict one copy on X
            ay_te <- rbind(ay_te, tmp[V15==fs])
        }  else{}
    }
}

a_out <- rbind(a_auto_te, ax_te) %>% rbind(.,ay_te)

total_copy <- a_out[te_source %like% 'hom',.N] * 2 + a_out[te_source %like% 'het',.N]

cat('========================\n\n')
cat('The total copy of TE in parents:', total_copy, '\n' )
cat('========================\n\n')
fwrite(a_out, outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep='\t')