#!/usr/bin/env Rscript
library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="Teflon output gt file")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--outfile", help="ouput file with alt frequency")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")
args <- parser$parse_args()
infile <- args$infile # infile='/data/home/ywang120/myData/PeterData/main_pipeline/idata/TE/dmel/fam1/genotypes/dmel_fam1_te_gt.txt'
pedfile <- args$pedfile # pedfile= '/data/home/ywang120/myData/PeterData/main_pipeline/idata/ped/dmel_fam1.ped'
outfile <- args$outfile
chrfile <- args$chrfile
###########################################

#MIN_SUPPORTING_READ = 3
MIN_FREQ = 0.2 # 97% area under MIN_TOTAL_READ=10
MIN_TOTAL_READ = 10 # presence reads + absence reads for each sample
NON_FOCAL_MAX_READ  = 1 # should < MIN_TOTAL_READ * MIN_FREQ
IMPURITY_SAMPLE_MAX = 3

################################################

chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]
Y_CHR <- chr[V1=='Y', V2]

###########################################

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
###########################################

aa <- fread(infile, header=F, sep="\t")
a_auto <- aa[V1 %in% AUTO_CHR]
#ambiguous_read_prop <- aa[,sum(.SD[,V12])/(sum(.SD[,V10+V11+V12])),by=V15]
idf <- a_auto[,c(1:4,14)] %>% unique


a_auto_cand <- data.table()

for(i in 1:nrow(idf)){
    print(i)
    chr = idf[i, V1]
    p5 = idf[i, V2]
    p3 = idf[i, V3]
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
    # all samples should meet minimum read depth
    if(!all(tmp[,V10+V11] >= MIN_TOTAL_READ)) next
    # all samples should meet maximum read depth, flaged as -9
    if(any(tmp[,V13]==-9)) next
    
    # Parental to be zero absence read
    if(!((tmp[V15==fs, V11] == 0)  & (tmp[V15==ms, V11]==0))) next 
    tmp[,absence_freq:=V11/(V10+V11)]
    # Allelic balance
    ab = max(tmp$absence_freq)
    if(tmp[absence_freq==ab,.N]>1) next # multiple max freqs
    if(!ab %between% c(MIN_FREQ, 1 - MIN_FREQ)) next
    # all non-focal samples should meet NON_FOCAL_MAX_READ
    if(!(tmp[V11 <= NON_FOCAL_MAX_READ, .N] == (sample_N-1))) next
    # number of impurity samples
    if(!(tmp[V11!=0,.N] <= (IMPURITY_SAMPLE_MAX + 1))) next
    a_auto_cand <- rbind(a_auto_cand, tmp)
}



#############################################


ax <- aa[V1 %in% X_CHR]
#ambiguous_read_prop <- aa[,sum(.SD[,V12])/(sum(.SD[,V10+V11+V12])),by=V15]

idf <- ax[,c(1:4,14)] %>% unique

a_x_cand <- data.table()
for(i in 1:nrow(idf)){
    print(i)
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
    # all females should meet minimum read depth on X
    if(!all(tmp[V15 %in% female_samples,V10+V11] >= MIN_TOTAL_READ)) next
    # all males should meet minimum read depth/2 on X
    if(!all(tmp[V15 %in% male_samples,V10+V11] >= MIN_TOTAL_READ/2)) next
    # all samples should meet maximum read depth, flaged as -9
    if(any(tmp[,V13]==-9)) next
    
    # mother and father to be zero supporting read
    if(!((tmp[V15==fs, V11] == 0)  & (tmp[V15==ms, V11]==0))) next
    tmp[,absence_freq:=V11/(V10+V11)]
    # Allelic balance
    ab = max(tmp$absence_freq)
    if(tmp[absence_freq==ab,.N]>1) next # multiple max freqs
    j <- which(tmp$absence_freq==ab)
    t_sample <- tmp[j,V15]
    if(t_sample %in% female_samples){
        if(!ab %between% c(MIN_FREQ, 1 - MIN_FREQ)) next
    }else {
       if(ab!=1) next # frequency need to be 1 in males
    }
    # all non-focal samples should meet NON_FOCAL_MAX_READ
    if(!(tmp[V11 <= NON_FOCAL_MAX_READ, .N] == (sample_N -1))) next
    # number of impurity samples
    if(!(tmp[V11!=0,.N] <= (IMPURITY_SAMPLE_MAX + 1))) next
    a_x_cand <- rbind(a_x_cand, tmp)
}

#############################################################

ay <- aa[V1 %in% Y_CHR]
a_y_cand <- data.table()


#ambiguous_read_prop <- aa[,sum(.SD[,V12])/(sum(.SD[,V10+V11+V12])),by=V15]
if(ay[,.N]!=0){
    idf <- ay[,c(1:4,14)] %>% unique
    for(i in 1:nrow(idf)){
        print(i)
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
        # all males should meet minimum read depth/2 on X
        if(!all(tmp[V15 %in% male_samples,V10+V11] >= MIN_TOTAL_READ/2)) next
        # all male samples should meet maximum read depth, flaged as -9
        if(any(tmp[V15%in% male_samples, V13]==-9)) next
        
        # mother and father to be zero supporting read
        if(!((tmp[V15==fs, V11] == 0)  & (tmp[V15==ms, V11]==0))) next
        tmp[,absence_freq:=V11/(V10+V11)]
        # Allelic balance
        ab = max(tmp$absence_freq, na.rm=TRUE)
        if(tmp[absence_freq==ab,.N]>1) next # multiple max freqs
        j <- which(tmp$absence_freq==ab)
        t_sample <- tmp[j,V15]
        if(t_sample %in% female_samples){
            next
        }else {
            if(ab!=1) next # frequency need to be 1 in males
        }
        # all non-focal samples should meet NON_FOCAL_MAX_READ
        if(!(tmp[V11 <= NON_FOCAL_MAX_READ, .N] == (sample_N -1))) next
        # number of impurity samples
        if(!(tmp[V11!=0,.N] <= (IMPURITY_SAMPLE_MAX + 1))) next
        a_y_cand <- rbind(a_y_cand, tmp)
    }
}



a_out <- rbind(a_auto_cand, a_x_cand) %>% rbind(.,a_y_cand)

#######################################


######################
##  cross family check
# we don't distingush the insertion or deletion
fam = stringr::str_split(infile, '_')[[1]][2]
oo_list <- list.files(pattern = '^M..genotypes.txt$', path='../..', recursive=TRUE, full.name = TRUE)

oo_list <- oo_list[!stringr::str_detect(oo_list, paste0('fam',fam))]

cat('The following files will be used for cross-family check:\n')
cat(paste(oo_list, collapse='\n'))
cat('\n')

if(a_out[,.N]!=0){
    idf <- a_out[,c(1:4,14)] %>% unique
    idf$presence_in_other_fam <- 0
    for(i in 1:nrow(idf)){
        chr = idf[i,V1]
        p5 = idf[i, V2]
        p3 = idf[i,V3]
        tid = idf[i,V4]
        tuid = idf[i, V14]
        for(j in oo_list){
            oo <- fread(j, header = FALSE, sep='\t')
            tmp <- oo[V1==chr & V2==p5 &  V3==p3 &V4==tid & V14==tuid]
            if(tmp[,.N]==0){next}
            else{
                idf[i,presence_in_other_fam:=1]
                cat('The following candidate presents in sample: ', j,'\n')
                cat(paste(c(chr, p5,p3,tid,tuid), collapse='\t'), '\n')
                break
            }
        }
    }
    a_out <- left_join(a_out, idf, by = c('V1','V2','V3','V4','V14'))
    a_out <- a_out[presence_in_other_fam==0]
}

fwrite(a_out, outfile, col.names=TRUE, row.names=FALSE, sep='\t', quote=FALSE)

