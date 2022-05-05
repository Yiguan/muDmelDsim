#!/usr/bin/env Rscript

library("data.table")
library("dplyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="input file with relevant annotation ")
parser$add_argument("--outfile", help="ouput file with condidate donovo")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="chromosome to be analysed")
parser$add_argument("--mindp", help="minimum read depth", default=10, type='integer')
parser$add_argument("--maxdp", help="maximum read depth", default=150, type='integer')
# parser$add_argument("--mingq", help="minimum genotype quality", default=50, type='integer')
parser$add_argument("--minab", help="allelic balance", default=0.15, type='double')
parser$add_argument("--minalt", help="minimum read count for the alternates", default=5, type='integer')

argv <- parser$parse_args()

MIN_DP = argv$mindp 
MAX_DP = argv$maxdp
##MIN_GQ = argv$mingq
AB     = argv$minab
MIN_ALT = argv$minalt
infile = argv$infile
outfile = argv$outfile
pedfile = argv$pedfile
chrfile = argv$chrfile
# MIN_DP = 10
# MAX_DP = 150
# MIN_GQ = 50 
# MIN_ALT = 5
# AB = 0.15

#########################################

chr <- fread(chrfile, header=F, sep=',')
auto_chr <- chr[V1=='AUTO', V2]
X_chr = chr[V1=='X', V2]
Y_chr = chr[V1=='Y', V2]
X_male_prop = 0.5 # denovo in male offspring should have >=100% ALT reads,
# but that's too strict for indels. Some may regard as soft-clip
# we choose 0.5, which will remove the candidate of  the same proportion as that in autosomes


#########################################
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
sample_n <- nrow(ped) + 2
#########################################

aa <- fread(infile, header=T, sep="\t")
names(aa) <- stringr::str_replace_all(names(aa), ms, "mother")
names(aa) <- stringr::str_replace_all(names(aa), fs, "father")

# autosome filtering
a_auto <- aa[CHROM %in% auto_chr]
# filter parents:  depth, genotype quality, purity
a_auto <- a_auto[mother.NR %between% c(MIN_DP, MAX_DP)  & mother.NV==0 &
                 father.NR %between% c(MIN_DP, MAX_DP)  & father.NV==0]

# filter non-focal unique sample
a_auto$zeroAlt <- a_auto[,rowSums(.SD==0),.SDcols=grep(".*NV$", names(a_auto), value=TRUE)]
a_auto <- a_auto[zeroAlt==(sample_n - 1)] # do not allow impurity samples

# filter alt count
a_auto$altCount<- a_auto[,rowSums(.SD),.SDcols=grep(".*NV$", names(a_auto), value=TRUE)]
a_auto <- a_auto[altCount >= MIN_ALT]

# filter allelic balance
alt1 <- a_auto[,.SD,.SDcols=grep(".*NV$", names(a_auto), value=TRUE)] %>% as.matrix()
dp1 <-  a_auto[,.SD,.SDcols=grep(".*NR$", names(a_auto), value=TRUE)] %>% as.matrix()
ab1 <- alt1/dp1
denovo_auto <- a_auto[apply(ab1,1,max) %between% c(AB, 1-AB),]

#######################################################################################
# X chromosome
ax <- aa[CHROM %in% X_chr]

ax <- ax[mother.NR %between% c(MIN_DP, MAX_DP)  & mother.NV==0 &
         father.NV == 0]

mv_male <- ax
mv_female <- ax[father.NR %between% c(MIN_DP/2, MAX_DP/2) ]



# For male offspring
mv_male$zeroAlt <- mv_male[,rowSums(.SD==0),.SDcols=grep(".*NV$", names(mv_male), value=TRUE)]
mv_male <- mv_male[zeroAlt==(sample_n - 1)]
mv_male$altCount <- mv_male[,rowSums(.SD),.SDcols=grep(".*NV$", names(mv_male), value=TRUE)]
mv_male <- mv_male[altCount >= MIN_ALT ]

denovo_X_male <- data.table()
if(nrow(mv_male)!=0){
    for(i in 1:nrow(mv_male))
    {
        tmp <- mv_male[i,]
        nr <- tmp[,.SD, .SDcols=c("CHROM","POS","REF","ALT", grep(".*NR",names(tmp), value=TRUE))] %>% 
            tidyr::gather("sid","NR",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NR")
        nv <- tmp[,.SD, .SDcols=c("CHROM","POS", "REF","ALT",grep(".*NV",names(tmp), value=TRUE))] %>% 
            tidyr::gather("sid","NV",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NV")
        rv <- nr %>% dplyr::left_join(.,nv[,c("sid","NV")], by = "sid") %>% mutate(ab=NV/NR) %>%
            filter(ab!=0)
        if(nrow(rv)!=1) stop("Impurity samples!")
        if(rv[1,"sid"] %in% female_offspring) next
        if(rv[1,"ab"] < X_male_prop) next 
        denovo_X_male <- rbind(denovo_X_male, tmp)
    }
}
# For female offspring
mv_female$zeroAlt <- mv_female[,rowSums(.SD==0),.SDcols=grep(".*NV$", names(mv_female), value=TRUE)]
mv_female <- mv_female[zeroAlt==13]
mv_female$altCount <- mv_female[,rowSums(.SD),.SDcols=grep(".*NV$", names(mv_female), value=TRUE)]
mv_female <- mv_female[altCount >= MIN_ALT ]

denovo_X_female <- data.table()
if(nrow(mv_female) !=0){
    for(i in 1:nrow(mv_female))
    {
        tmp <- mv_female[i,]
        nr <- tmp[,.SD, .SDcols=c("CHROM","POS", "REF","ALT", grep(".*NR",names(tmp), value=TRUE))] %>% 
            tidyr::gather("sid","NR",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NR")
        nv <- tmp[,.SD, .SDcols=c("CHROM","POS", "REF","ALT", grep(".*NV",names(tmp), value=TRUE))] %>% 
            tidyr::gather("sid","NV",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NV")
        rv <- nr %>% dplyr::left_join(.,nv[,c("sid","NV")], by = "sid") %>% mutate(ab=NV/NR) %>%
            filter(ab!=0)
        if(nrow(rv)!=1) stop("Impurity samples!")
        if(rv[1,"sid"] %in% male_offspring) next
        if(rv[1,"ab"] < AB || rv[1,"ab"]>1-AB) next 
        denovo_X_female <- rbind(denovo_X_female, tmp)
    }
}
denovo_X <- rbind(denovo_X_male, denovo_X_female)




#######################################################################################
# Y chromosome
ay <- aa[CHROM %in% Y_chr]

ay <- ay[father.NR %between% c(MIN_DP/2, MAX_DP/2)  & mother.NV==0 &
         father.NV == 0]
denovo_ay <- data.table()
# For male offspring
if(nrow(ay) !=0){
    ay$zeroAlt <- ay[,rowSums(.SD==0),.SDcols=grep(".*NV$", names(ay), value=TRUE)]
    ay <- ay[zeroAlt==(sample_n - 1)]
    ay$altCount <- ay[,rowSums(.SD),.SDcols=grep(".*NV$", names(ay), value=TRUE)]
    ay <- ay[altCount >= MIN_ALT ]
    if(nrow(ay)!=0){
        for(i in 1:nrow(ay))
        {
            tmp <- ay[i,]
            nr <- tmp[,.SD, .SDcols=c("CHROM","POS", "REF","ALT", grep(".*NR",names(tmp), value=TRUE))] %>% 
                tidyr::gather("sid","NR",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NR")
            nv <- tmp[,.SD, .SDcols=c("CHROM","POS", "REF","ALT", grep(".*NV",names(tmp), value=TRUE))] %>% 
                tidyr::gather("sid","NV",5:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NV")
            rv <- nr %>% dplyr::left_join(.,nv[,c("sid","NV")], by = "sid") %>% mutate(ab=NV/NR) %>%
                filter(ab!=0)
            if(nrow(rv)> 1) stop("Impurity samples!")
            if(rv[1,"sid"] %in% female_offspring) next
            if(rv[1,"ab"] < X_male_prop) next 
            denovo_ay <- rbind(denovo_ay, tmp)
        }
    }
}

denovo_XY <- rbind(denovo_X, denovo_ay)


denovo <- rbind(denovo_auto, denovo_XY)
denovo$zeroAlt <- NULL
denovo$altCount <- NULL


# ## format output
# a12 <- denovo[,.SD,.SDcols=c("CHROM","POS", grep(".*NV",names(denovo), value=TRUE))] %>%
#     tidyr::gather("sid","NV",3:ncol(.)) %>% tidyr::extract("sid","sid","(.*).NV") %>%
#     filter(NV!=0)
# outdt <- data.table()
# for(j in 1:nrow(a12))
# {
#     chr <- a12[j,"CHROM"]
#     pos <- a12[j,"POS"]
#     cmd1 = stringr::str_interp("awk '{if($1==\"${chr}\" && $2==${pos}){print $4}}' allsample_indel_biallelic.recode.vcf")
#     ref <- system(cmd1,intern = TRUE)
#     cmd2 = stringr::str_interp("awk '{if($1==\"${chr}\" && $2==${pos}){print $5}}' allsample_indel_biallelic.recode.vcf")
#     alt <- system(cmd2,intern = TRUE)
#     tmp <- data.table("CHROM"=chr, "POS"=pos, "REF"=ref, "ALT"=alt, "sid"=a12[j,"sid"])
#     outdt <- rbind(outdt, tmp)
# }

fwrite(denovo, outfile, row.names=F, col.names=T, sep = "\t", quote=F)