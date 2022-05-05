#!/usr/bin/env Rscript

library("data.table")
library("dplyr")
library("tidyr")
library("ggplot2")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="input file with phased vcf-liked body, seperated by comma and with header")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="comma-separated file about chromosomes")

args <- parser$parse_args()
infile <- args$infile # infile = 'dmel_fam1.snp.hardFiltered.biallelic.GQ20.phased.gz'
pedfile <- args$pedfile
chrfile <- args$chrfile


################################################

chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]

CHR <- c(AUTO_CHR, X_CHR)
#CHR = c('2L','2R','3L','3R', 'X') # TO DO: what if reference genome has different names!
#####################################################

## fix 11111222211111 -> 111111111111
## 111221122 -> 33333333
# fix based on flanking regions;
# if flanking regions too short, treat as missing
fixBase <- function(ivec){
    MAX_MISS = 4
    FLANK_MATCH = 4 # on each side
    #ivec <- m1$father_h1
    t1 <- rle(ivec)
    tt <- data.table("len"=t1$lengths, "vv"=t1$values)
    for(r in which(tt$len<=MAX_MISS)){
        pre <- r-1
        nex <- r+1
        if(pre<1 || nex>nrow(tt)){
            tt$vv[r] <- 3
        }
        else if(tt$len[pre]>=FLANK_MATCH && tt$len[nex]>=FLANK_MATCH){
            tt$vv[r] <- tt$vv[r]%%2 + 1 # 1 -> 2 or 2 -> 1
        }else{
            tt$vv[r] <- 3 # regard missing
        }
    }
    re <- rep(tt$vv, tt$len)
    return(re)
}

# hap plot by SNPs
plotSNP <- function(plotdf, chr="2L", legend_title="maternal"){
    plotdt <- plotdf[CHROM==chr]
    plotdt$snpid <- 1:nrow(plotdt)
    plotdt <- plotdt %>% gather("ss","hh",grep("_source",names(plotdt), value=TRUE))
    plotdt$ss <- factor(plotdt$ss)
    ggplot(data=plotdt) + geom_rect(aes(xmin=snpid-0.5,xmax=snpid+0.5,ymin=as.integer(ss)-0.5, 
                                        ymax=as.integer(ss)+0.5, fill=as.factor(hh), color=as.factor(hh))) +
        scale_y_continuous(breaks=as.integer(unique(plotdt$ss)), expand=c(0,0),
                           labels=stringr::str_remove(levels(unique(plotdt$ss)),"_source")) +
        scale_x_continuous(name='SNP markers', expand=c(0,0)) + theme_bw() +
        scale_fill_discrete(name=legend_title, labels=c("hap1","hap2")) +
        scale_color_discrete(name=legend_title, labels=c("hap1","hap2"))
}
###
getHap_breakpoint <- function(plotdf, chr = "2L", parents = "mother"){
    plotdt = plotdf[CHROM==chr]
    ret_list <- vector(mode = "list",length=2)
    bpdf <- data.table()
    hapdf <- data.table()
    for(s in grep("_source",names(plotdt), value = TRUE)){
        print(s)
        t1 <- rle(plotdt[[s]])
        t2 <- data.table(ll=t1$lengths,vv=t1$values) %>% mutate(s_end=cumsum(ll)) %>%
            mutate(s_start=s_end-ll+1)
        t2$bp_start <- plotdt$POS[t2$s_start]
        t2$bp_end <- plotdt$POS[t2$s_end]
        t2$bp_len <- t2$bp_end - t2$bp_start + 1
        t2$parent <- parents
        t2$sid <- stringr::str_remove(s,"_source")
        t2$chr <- chr
        t2 <- t2[,c(9,10,8,2,4,3,1,5,6,7)]
        names(t2) <- c("sample","chr","parent","hapid","snp_s","snp_e","snp_n","bp_s","bp_e","bp_n")
        hapdf <- rbind(hapdf,t2)
        ## breakpoint
        # remove NA and compact continous hapblocks
        noNA_hap <- t2[hapid!=3]
        m1 <- rle(noNA_hap[,hapid])
        m2 <- data.table(ll=m1$lengths,vv=m1$values) %>% mutate(lineid=cumsum(ll))
        m_fix <- m2[ll!=1]
        if(m_fix[,.N]!=0){
            for(i in 1:nrow(m_fix)){
                fix_line_start  <- m_fix[i,lineid] - m_fix[i,ll] + 1
                fix_line_end <- m_fix[i,lineid]
                fix_snp_s <- noNA_hap[fix_line_start,snp_s]
                fix_snp_e <- noNA_hap[fix_line_end, snp_e]
                fix_snp_n <- sum(noNA_hap[fix_line_start:fix_line_end,snp_n])
                fix_bp_s <- noNA_hap[fix_line_start,bp_s]
                fix_bp_e <- noNA_hap[fix_line_end, bp_e]
                fix_bp_n <- sum(noNA_hap[fix_line_start:fix_line_end,bp_n])
                # make duplicate lines
                noNA_hap[fix_line_start:fix_line_end,snp_s:=fix_snp_s]
                noNA_hap[fix_line_start:fix_line_end,snp_e:=fix_snp_e]
                noNA_hap[fix_line_start:fix_line_end,snp_n:=fix_snp_n]
                noNA_hap[fix_line_start:fix_line_end,bp_s:=fix_bp_s]
                noNA_hap[fix_line_start:fix_line_end,bp_e:=fix_bp_e]
                noNA_hap[fix_line_start:fix_line_end,bp_n:=fix_bp_n]
            }
        }
        noNA_hap <- unique(noNA_hap) # if show hap ignoring NA, use this datatable instead of 'hapdf'
        if(noNA_hap[,.N]>1){ # in case no recombination
            breakpoint <- data.table("sample" = stringr::str_remove(s,"_source"),
                                     "chr" = unique(plotdt$CHROM),
                                     "parent" = parents,
                                     "breakpoint_s" = noNA_hap[1:(.N-1),bp_e],
                                     "breakpoint_e" = noNA_hap[2:.N, bp_s],
                                     "snp_s" = noNA_hap[1:(.N-1),snp_e],
                                     "snp_e" = noNA_hap[2:.N,snp_s])
            bpdf <- rbind(bpdf, breakpoint)
        }
    }
    ret_list[[1]] <- hapdf
    ret_list[[2]]<- bpdf
    return(ret_list)
}   
###
plotHap <- function(hapdf, legend_title = "maternal"){
    hapdf$sample <- as.factor(hapdf$sample)
    ggplot(hapdf) + geom_rect(aes(xmin=bp_s,xmax=bp_e,ymin=as.integer(hapdf$sample)-0.5,
                                 ymax=as.integer(hapdf$sample)+0.5, fill=as.factor(hapid)), color=NA) +
        scale_y_continuous(breaks=as.integer(unique(hapdf$sample)), expand=c(0,0),
                           labels=levels(unique(hapdf$sample))) +
        scale_x_continuous(name='genomic location (bp)',expand=c(0,0)) + theme_bw() +
        scale_fill_discrete(name=legend_title, labels=c("hap1","hap2"))
}

###################################################################
ped <- fread(pedfile, header=F, sep='\t')
ms <- unique(ped$V4)
if(length(ms)!=1){
    stop('PED: samples are NOT from the same mother!')
}
fs <- unique(ped$V3)
if(length(fs)!=1){
    stop('PED: samples are NOT from the same father!')
}
#########################################


aa <- fread(infile, header=T, sep=",")
# for each site, all samples should be genotyped and phased (homo is considered phased)
all_phased <- which(apply(aa,1,function(x) length(grep("/",x)))==0) 
aa <- aa[all_phased]
aa <- aa %>% gather("ss","gt",names(aa)[5:ncol(aa)]) %>% extract("gt",c("h1","h2"), "(.).(.)") %>% as.data.table

aa$h1 <- as.integer(aa$h1)
aa$h2 <- as.integer(aa$h2)
aa <- aa %>% gather("h12","gcode",c("h1","h2")) %>%
  unite("s_hap", ss:h12, remove=T) %>% spread(s_hap, gcode) %>% as.data.table

# rename parental samples
new_names <- stringr::str_replace_all(names(aa), paste0(ms,'_h'), 'mother_hap')
new_names <- stringr::str_replace_all(new_names, paste0(fs,'_h'), 'father_hap')
names(aa) <- new_names

########################################
out_haps <- data.table()
out_breakpoints <- data.table()

##############
## maternal haps
am <- aa[,.SD,.SDcols=c("CHROM","POS","ALT","REF", "mother_hap1","mother_hap2",
                        grep("_h1",names(aa), value=TRUE))]  
# filter het site in parents
am <- am[mother_hap1!=mother_hap2]
for(ss in grep("_h1", names(am), value=TRUE)){
    s_check <- stringr::str_replace(ss,"_h1","_source")
    am[,eval(s_check):=2] # from second parental hap
    am[get(ss)==mother_hap1,eval(s_check):=1] # from first parenta l hap
}
am <- am[,.SD, .SDcol=c("CHROM","POS",grep(".*_source$",names(am), value=TRUE))]

# fix base
for(chr in CHR){
    am[CHROM==chr, 3:ncol(am)] <- sapply(am[CHROM==chr,3:ncol(am)], fixBase) %>% as.data.table
}

# plot haplotypes using SNP markers
for(chr in CHR){
    pp <- plotSNP(am, chr = chr, legend_title = "maternal")
    ggsave(paste0(chr,'_maternal_snp.pdf'), pp, width=11,
           height=6, units="in", dpi=100)
    mo_test <- getHap_breakpoint(am, chr = chr, parents = "mother")
    out_haps <- rbind(out_haps, mo_test[[1]])
    out_breakpoints <- rbind(out_breakpoints, mo_test[[2]])
    ww <- plotHap(mo_test[[1]], legend_title = "maternal")
    ggsave(paste0(chr,'_maternal_hap.pdf'), ww, width=11,
        height=6, units="in", dpi=100)
}

##############
## paternal haps
af <- aa[,.SD,.SDcols=c("CHROM","POS","ALT","REF", "father_hap1","father_hap2",
                        grep("_h2",names(aa), value=TRUE))]  
# filter het site in parents
af <- af[father_hap1!=father_hap2]
for(ss in grep("_h2", names(af), value=TRUE)){
    s_check <- stringr::str_replace(ss,"_h2","_source")
    af[,eval(s_check):=2] # from second parental hap
    af[get(ss)==father_hap1,eval(s_check):=1] # from first parenta l hap
}
af <- af[,.SD, .SDcol=c("CHROM","POS",grep(".*_source$",names(af), value=TRUE))]

# fix base
for(chr in CHR){
    af[CHROM==chr, 3:ncol(af)] <- sapply(af[CHROM==chr,3:ncol(af)], fixBase) %>% as.data.table
}

# plot haplotypes using SNP markers
for(chr in CHR){
    pp <- plotSNP(af, chr = chr, legend_title = "paternal")
    ggsave(paste0(chr,'_paternal_snp.pdf'), pp, width=11,
           height=6, units="in", dpi=100)
    fa_test <- getHap_breakpoint(af, chr = chr, parents = "father")
    out_haps <- rbind(out_haps, fa_test[[1]])
    out_breakpoints <- rbind(out_breakpoints, fa_test[[2]])
    ww <- plotHap(fa_test[[1]], legend_title = "paternal")
    ggsave(paste0(chr,'_paternal_hap.pdf'), ww, width=11,
        height=6, units="in", dpi=100)
}

## Output haps and breakpoints

fwrite(out_haps, "hap_range.txt", col.names=T, row.names=F, sep='\t', quote=F)
fwrite(out_breakpoints, "breakpoints_range.txt", col.names=T, row.names=F, sep='\t', quote=F)




