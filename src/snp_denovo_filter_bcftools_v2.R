#!/usr/bin/env Rscript
library("data.table")
library("dplyr")
library("tidyr")
library("argparse")

## version log
# based on 'callable_by_simulation_bcftools.R'
# consider variants in other families when estimating the callable sites



parser <- ArgumentParser()
parser$add_argument("--infile", help="input file with relevant annotation table.gz")
parser$add_argument("--outfile", help="ouput file with condidate donovo")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="comma separated file, telling whether auto or X chromosome")
parser$add_argument("--famdir", help="eg: /data/home/ywang120/myData/bbsrc/idata/denovo_mutations/dmel")
parser$add_argument("--mindp", help="minimum read depth", default=10, type='integer')
parser$add_argument("--maxdp", help="maximum read depth", default=150, type='integer')
parser$add_argument("--mingq", help="minimum genotype quality", default=50, type='integer')
parser$add_argument("--minab", help="allelic balance", default=0.2, type='double')
parser$add_argument("--minalt", help="minimum read count for the alternates", default=5, type='integer')

argv <- parser$parse_args()

MIN_DP = argv$mindp 
MAX_DP = argv$maxdp
MIN_GQ = argv$mingq
AB     = argv$minab
MIN_ALT = argv$minalt
infile = argv$infile
outfile = argv$outfile
pedfile = argv$pedfile
chrfile = argv$chrfile
VARIANTS_FAM_DIR = argv$famdir
#################################################
# MIN_DP = 10
# MAX_DP = 150
# MIN_GQ = 50
# AB = 0.20
# MIN_ALT = 5
# when filter_impurity = FALSE, it means filter non-focal samples just based on GT
filter_impurity = TRUE

# When the non-focal samples with AD alt >= ad_min_for_impurity,
# the sample is then considered as 'impurity sample'
ad_min_for_impurity = 1 #integer [1, Inf]

# When the number of impurity samples >= this value,
# the loci will not be considered as a deonvo
max_impurity_samples = 2 # integer [0,Inf]

# AUTO_CHR = c('2L','2R','3L','3R')
# X_CHR = c('X')

SPECIES = stringr::str_split(basename(pedfile),'_')[[1]][1]
FAM = stringr::str_split(stringr::str_split(basename(pedfile),'_')[[1]][2],'\\.')[[1]][1]



####################################################
# check AD
checkAD <- function(ad){ #ad='20,0' ; ad='0,20,17', ad='5,0,0'
    tt <- stringr::str_split(ad,',')
    alt_c <- tt[2:length(tt)]
    ret_v <- lapply(tt, function(x) x[1]!='0' & all(x[2:length(x)]=='0')) %>% unlist
    return(ret_v)
}
#####################################################

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
SAMPLE_N <- nrow(ped)+2


chr <- fread(chrfile, header=F, sep=',')
AUTO_CHR <- chr[V1=='AUTO', V2]
X_CHR <- chr[V1=='X', V2]
Y_CHR <- chr[V1=='Y', V2]
#####################################################

aa <- fread(infile, header = TRUE, sep = "\t")
names(aa) <- stringr::str_replace_all(names(aa), ms, "mother")
names(aa) <- stringr::str_replace_all(names(aa), fs, "father")


###################### autosome ##############################
## Filter parents
# autosome filtering
a_auto <- aa[CHROM %in% AUTO_CHR]

a_auto <- a_auto[(!mother.GT %like% "\\*") & (!father.GT %like% "\\*")] # parental should no contain deletions
# filter parents: Read depth and genotype quality
a_auto <- a_auto[mother.DP %between% c(MIN_DP, MAX_DP) & mother.GQ >=MIN_GQ & 
                father.DP %between% c(MIN_DP, MAX_DP) & father.GQ >= MIN_GQ]
# parents should be purity sample for reference allele. AD=x:0
mv_auto <- a_auto[checkAD(mother.AD) & checkAD(father.AD),]


mvdt_auto <- data.table()
for(i in 1:nrow(mv_auto)){
    print(i)
    tmp <- mv_auto[i,]
    ref <- tmp[1,REF]
    refGT <- paste0(ref,"/",ref)
    alt <- stringr::str_split(tmp[1,ALT],",")[[1]]
    # only consider het offspring
    altGT <- paste0(ref,"/",alt) # assume and checked all GT, ref came as first
    
    tmp.gt <- tmp[,.SD,.SDcols=c("CHROM","POS","REF","ALT", grep(".*GT$", names(mv_auto), value=TRUE))]
    tmp.gt <- tmp.gt %>% gather("sid","igt",5:(SAMPLE_N+4)) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.dp <- tmp[,.SD,.SDcols= grep(".*DP$", names(mv_auto), value=TRUE)]
    tmp.dp <- tmp.dp %>% gather("sid","idp",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.gq <- tmp[,.SD,.SDcols= grep(".*GQ$", names(mv_auto), value=TRUE)]
    tmp.gq <- tmp.gq %>% gather("sid","igq",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.ad <- tmp[,.SD,.SDcols= grep(".*AD$", names(mv_auto), value=TRUE)]
    tmp.adf <- tmp[,.SD,.SDcols= grep(".*ADF$", names(mv_auto), value=TRUE)]
    tmp.adr <- tmp[,.SD,.SDcols= grep(".*ADR$", names(mv_auto), value=TRUE)]
    if(length(alt)==1){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1"),"(.*),(.*)")
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1"),"(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1"),"(.*),(.*)") 
    }
    if(length(alt)==2){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1","ad_alt2"),"(.*),(.*),(.*)")
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1","adf_alt2"),"(.*),(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1","adr_alt2"),"(.*),(.*),(.*)") 
    }
    if(length(alt)>2){
        cat("> 2 alternative allele, skipped!\n")
        next
    }
    
    tmp_dt <- tmp.gt %>% left_join(.,tmp.gq, by = "sid") %>% left_join(.,tmp.dp, by="sid") %>%
        left_join(.,tmp.ad, by="sid") %>% left_join(.,tmp.adr, by="sid") %>% 
        left_join(.,tmp.adf, by="sid") 
    tmp_dt <- tmp_dt %>% mutate_at(names(tmp_dt)[7:ncol(tmp_dt)], as.integer) %>% 
        as.data.table  
    
    jalt = 1
    for(mu in altGT){
        adcol = paste0("ad_alt",jalt)
        adrcol = paste0("adr_alt",jalt)
        adfcol = paste0("adf_alt",jalt)
        jalt = jalt + 1
        cdt <- tmp_dt[igt==mu]
        #"Candidate mutation must not be present in any other unrelated sample"
        # non-unique 'mutation'; determined by GT, 
        if(cdt[,.N] !=1) next
        # filter by AD in other samples
        if(filter_impurity){
            pp <- tmp_dt[sid!=cdt$sid,adcol, with=F]
            names(pp) <- "ad_i"
            if(pp[ad_i > ad_min_for_impurity, .N]>0) next
            if(pp[ad_i > 0 & ad_i <= ad_min_for_impurity,.N] > max_impurity_samples) next
            
        }
        # focal child GQ and DP filter
        if(!((cdt[1,igq] >= MIN_GQ) && (cdt[1, idp] %between% c(MIN_DP, MAX_DP)))) next
        # focal child ALT count filter
        if(!(cdt[[1,adcol]] >= MIN_ALT)) next
        # "Candidate mutations must be present on both the forward and reverse strand in the offspring (ADF, ADR > 0)."
        if(!(cdt[[1,adrcol]] > 0 && cdt[[1,adfcol]] > 0)) next
        # "Candidate mutation must not have low allelic depth in the offspring" AB
        alti <- cdt[[1, adcol]]
        ab <- alti/(alti + cdt[1,ad_ref])
        if(!(ab %between% c(AB, 1-AB))) next
        mvdt_auto <- rbind(mvdt_auto, cdt[1,1:5])
    }
}


############################### X chromosomes ##################

# whatever male or female offspring, parents should be:
ax <- aa[CHROM %in% X_CHR]
ax <- ax[checkAD(mother.AD) & checkAD(father.AD)]

# asterisk indicates deletion in gatk genotype
ax <- ax[(!mother.GT %like% "\\*") & (!father.GT %like% "\\*")]

# male offspring only need filter mother
mv_x_male <- ax[mother.DP %between% c(MIN_DP, MAX_DP) & mother.GQ >=MIN_GQ]
# female offspring need to filter mother and father
mv_x_female <- ax[mother.DP %between% c(MIN_DP, MAX_DP) & mother.GQ >=MIN_GQ & 
              father.DP %between% c(MIN_DP/2, MAX_DP/2) & father.GQ >=MIN_GQ]

# MALE OFFSPRING
mvdt_male <- data.table()
for(i in 1:nrow(mv_x_male)){
    print(i)
    tmp <- mv_x_male[i,]
    ref <- tmp[1,REF]
    refGT <- paste0(ref,"/",ref)
    alt <- stringr::str_split(tmp[1,ALT],",")[[1]]
    #altGT <- paste0(alt,"/",alt) # for males, should be alt homo !!!!!!!
    
    tmp.gt <- tmp[,.SD,.SDcols=c("CHROM","POS","REF","ALT", grep(".*GT$", names(mv_x_male), value=TRUE))]
    tmp.gt <- tmp.gt %>% gather("sid","igt",5:(SAMPLE_N+4)) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.dp <- tmp[,.SD,.SDcols= grep(".*DP$", names(mv_x_male), value=TRUE)]
    tmp.dp <- tmp.dp %>% gather("sid","idp",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.gq <- tmp[,.SD,.SDcols= grep(".*GQ$", names(mv_x_male), value=TRUE)]
    tmp.gq <- tmp.gq %>% gather("sid","igq",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.ad <- tmp[,.SD,.SDcols= grep(".*AD$", names(mv_x_male), value=TRUE)]
    tmp.adf <- tmp[,.SD,.SDcols= grep(".*ADF$", names(mv_auto), value=TRUE)]
    tmp.adr <- tmp[,.SD,.SDcols= grep(".*ADR$", names(mv_auto), value=TRUE)]
    if(length(alt)==1){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1"),"(.*),(.*)") 
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1"),"(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1"),"(.*),(.*)") 
    }
    if(length(alt)==2){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1","ad_alt2"),"(.*),(.*),(.*)") 
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1","adf_alt2"),"(.*),(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1","adr_alt2"),"(.*),(.*),(.*)") 
    }
    if(length(alt)>2){
        cat("> 2 alternative allele, skipped!\n")
        next
    }
    
    tmp_dt <- tmp.gt %>% left_join(.,tmp.gq, by = "sid") %>% left_join(.,tmp.dp, by="sid") %>%
        left_join(.,tmp.ad, by="sid") %>% left_join(.,tmp.adr, by="sid") %>% 
        left_join(.,tmp.adf, by="sid") 
    tmp_dt <- tmp_dt %>% mutate_at(names(tmp_dt)[7:ncol(tmp_dt)], as.integer) %>% 
        as.data.table  
    
    jalt = 1
    for(mu in alt){
        adcol = paste0("ad_alt",jalt)
        adrcol = paste0("adr_alt",jalt)
        adfcol = paste0("adf_alt",jalt)
        jalt = jalt + 1
        cdt <- tmp_dt[igt==mu]
        #"Candidate mutation must not be present in any other unrelated sample"
        # non-unique 'mutation'; determined by GT, 
        if(cdt[,.N] !=1) next
        # if candidate is female, pass
        # if(cdt[1,sid] %in% female_offspring) next #
        # filter by AD in other samples
        if(filter_impurity){
            pp <- tmp_dt[sid!=cdt$sid,adcol, with=F]
            names(pp) <- "ad_i"
            if(pp[ad_i > ad_min_for_impurity, .N]>0) next
            if(pp[ad_i > 0 & ad_i <= ad_min_for_impurity,.N] > max_impurity_samples) next
            
        }
         # "Candidate mutations must be present on both the forward and reverse strand in the offspring (ADF, ADR > 0)."
        if(!(cdt[[1,adrcol]] > 0 && cdt[[1,adfcol]] > 0)) next
        # focal child GQ and DP filter
        if(!((cdt[1,igq] >= MIN_GQ) && (cdt[1, idp] %between% c(MIN_DP/2, MAX_DP/2)))) next
        # focal child ALT count filter
        if(!(cdt[[1,adcol]] >= MIN_ALT)) next # We don't half minimum alt even if in males
        # focal child pure sample
        alti <- cdt[[1, adcol]]
        ab <- alti/(alti + cdt[1,ad_ref])
        if(ab!=1) next
        mvdt_male <- rbind(mvdt_male, cdt[1,1:5])
    }
}

# FEMALE OFFSPRING
mvdt_female <- data.table()
for(i in 1:nrow(mv_x_female)){
    print(i)
    tmp <- mv_x_female[i,]
    ref <- tmp[1,REF]
    refGT <- paste0(ref,"/",ref)
    alt <- stringr::str_split(tmp[1,ALT],",")[[1]]
    # only consider het offspring
    altGT <- paste0(ref,"/",alt)
    
    tmp.gt <- tmp[,.SD,.SDcols=c("CHROM","POS","REF","ALT", grep(".*GT$", names(mv_x_female), value=TRUE))]
    tmp.gt <- tmp.gt %>% gather("sid","igt",5:(SAMPLE_N+4)) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.dp <- tmp[,.SD,.SDcols= grep(".*DP$", names(mv_x_female), value=TRUE)]
    tmp.dp <- tmp.dp %>% gather("sid","idp",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.gq <- tmp[,.SD,.SDcols= grep(".*GQ$", names(mv_x_female), value=TRUE)]
    tmp.gq <- tmp.gq %>% gather("sid","igq",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
    tmp.ad <- tmp[,.SD,.SDcols= grep(".*AD$", names(mv_x_female), value=TRUE)]
    tmp.adf <- tmp[,.SD,.SDcols= grep(".*ADF$", names(mv_auto), value=TRUE)]
    tmp.adr <- tmp[,.SD,.SDcols= grep(".*ADR$", names(mv_auto), value=TRUE)]
    if(length(alt)==1){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1"),"(.*),(.*)") 
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1"),"(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1"),"(.*),(.*)") 
    }
    if(length(alt)==2){
        tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iad",c("ad_ref","ad_alt1","ad_alt2"),"(.*),(.*),(.*)") 
        tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadf",c("adf_ref","adf_alt1","adf_alt2"),"(.*),(.*),(.*)") 
        tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
            extract("iadr",c("adr_ref","adr_alt1","adr_alt2"),"(.*),(.*),(.*)") 
    }
    if(length(alt)>2){
        cat("> 2 alternative allele, skipped!\n")
        next
    }
    
    tmp_dt <- tmp.gt %>% left_join(.,tmp.gq, by = "sid") %>% left_join(.,tmp.dp, by="sid") %>%
        left_join(.,tmp.ad, by="sid") %>% left_join(.,tmp.adr, by="sid") %>% 
        left_join(.,tmp.adf, by="sid") 
    tmp_dt <- tmp_dt %>% mutate_at(names(tmp_dt)[7:ncol(tmp_dt)], as.integer) %>% 
        as.data.table  
    
    jalt = 1
    for(mu in altGT){
        adcol = paste0("ad_alt",jalt)
        adrcol = paste0("adr_alt",jalt)
        adfcol = paste0("adf_alt",jalt)
        jalt = jalt + 1
        cdt <- tmp_dt[igt==mu]
        #"Candidate mutation must not be present in any other unrelated sample"
        # non-unique 'mutation'; determined by GT, 
        if(cdt[,.N] !=1) next
        # if candidate is male, pass
        # if(cdt[1,sid] %in% male_offspring) next # impossible in males based on GT
        # filter by AD in other samples
        if(filter_impurity){
            pp <- tmp_dt[sid!=cdt$sid,adcol, with=F]
            names(pp) <- "ad_i"
            if(pp[ad_i > ad_min_for_impurity, .N]>0) next
            if(pp[ad_i > 0 & ad_i <= ad_min_for_impurity,.N] > max_impurity_samples) next
            
        }
         # "Candidate mutations must be present on both the forward and reverse strand in the offspring (ADF, ADR > 0)."
        if(!(cdt[[1,adrcol]] > 0 && cdt[[1,adfcol]] > 0)) next
        # focal child GQ and DP filter
        if(!((cdt[1,igq] >= MIN_GQ) && (cdt[1, idp] %between% c(MIN_DP, MAX_DP)))) next
        # focal child ALT count filter
        if(!(cdt[[1,adcol]] >= MIN_ALT)) next
        # "Candidate mutation must not have low allelic depth in the offspring" AB
        alti <- cdt[[1, adcol]]
        ab <- alti/(alti + cdt[1,ad_ref])
        if(!(ab %between% c(AB, 1-AB))) next
        mvdt_female <- rbind(mvdt_female, cdt[1,1:5])
    }
}



############################### Y chromosomes ##################

# whatever male or female offspring, parents should be:
ay <- aa[CHROM %in% Y_CHR]
mvdt_y <- data.table()
if(nrow(ay)!=0){
    ay <- ay[checkAD(father.AD)]
    # asterisk indicates deletion in gatk genotype
    ay <- ay[(!father.GT %like% "\\*")]
    # male offspring only need filter mother
    mv_y <- ay[father.DP %between% c(MIN_DP/2, MAX_DP/2) & father.GQ >=MIN_GQ]
    # MALE OFFSPRING
    for(i in 1:nrow(mv_y)){
        print(i)
        tmp <- mv_y[i,]
        ref <- tmp[1,REF]
        #refGT <- paste0(ref,"/",ref)
        alt <- stringr::str_split(tmp[1,ALT],",")[[1]]
        #altGT <- paste0(alt,"/",alt) # for males, should be alt homo !!!!!!!
        tmp.gt <- tmp[,.SD,.SDcols=c("CHROM","POS","REF","ALT", grep(".*GT$", names(mv_y), value=TRUE))]
        tmp.gt <- tmp.gt %>% gather("sid","igt",5:(SAMPLE_N+4)) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
        tmp.dp <- tmp[,.SD,.SDcols= grep(".*DP$", names(mv_y), value=TRUE)]
        tmp.dp <- tmp.dp %>% gather("sid","idp",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
        tmp.gq <- tmp[,.SD,.SDcols= grep(".*GQ$", names(mv_y), value=TRUE)]
        tmp.gq <- tmp.gq %>% gather("sid","igq",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>% as.data.table
        tmp.ad <- tmp[,.SD,.SDcols= grep(".*AD$", names(mv_y), value=TRUE)]
        tmp.adf <- tmp[,.SD,.SDcols= grep(".*ADF$", names(mv_y), value=TRUE)]
        tmp.adr <- tmp[,.SD,.SDcols= grep(".*ADR$", names(mv_y), value=TRUE)]
        if(length(alt)==1){
            tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iad",c("ad_ref","ad_alt1"),"(.*),(.*)") 
            tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iadf",c("adf_ref","adf_alt1"),"(.*),(.*)") 
            tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iadr",c("adr_ref","adr_alt1"),"(.*),(.*)") 
        }
        if(length(alt)==2){
            tmp.ad <- tmp.ad %>% gather("sid","iad",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iad",c("ad_ref","ad_alt1","ad_alt2"),"(.*),(.*),(.*)") 
            tmp.adf <- tmp.adf %>% gather("sid","iadf",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iadf",c("adf_ref","adf_alt1","adf_alt2"),"(.*),(.*),(.*)") 
            tmp.adr <- tmp.adr %>% gather("sid","iadr",1:SAMPLE_N) %>% extract("sid","sid","(.*)\\.") %>%
                extract("iadr",c("adr_ref","adr_alt1","adr_alt2"),"(.*),(.*),(.*)") 
        }
        if(length(alt)>2){
            cat("> 2 alternative allele, skipped!\n")
            next
        }
        
        tmp_dt <- tmp.gt %>% left_join(.,tmp.gq, by = "sid") %>% left_join(.,tmp.dp, by="sid") %>%
            left_join(.,tmp.ad, by="sid") %>% left_join(.,tmp.adr, by="sid") %>% 
            left_join(.,tmp.adf, by="sid") 
        tmp_dt <- tmp_dt %>% mutate_at(names(tmp_dt)[7:ncol(tmp_dt)], as.integer) %>% 
            as.data.table  
        
        jalt = 1
        for(mu in alt){
            adcol = paste0("ad_alt",jalt)
            adrcol = paste0("adr_alt",jalt)
            adfcol = paste0("adf_alt",jalt)
            jalt = jalt + 1
            cdt <- tmp_dt[igt==mu]
            #"Candidate mutation must not be present in any other unrelated sample"
            # non-unique 'mutation'; determined by GT, 
            if(cdt[,.N] !=1) next
            # if candidate is female, pass
            # if(cdt[1,sid] %in% female_offspring) next #
            # filter by AD in other samples
            if(filter_impurity){
                pp <- tmp_dt[sid!=cdt$sid,adcol, with=F]
                names(pp) <- "ad_i"
                if(pp[ad_i > ad_min_for_impurity, .N]>0) next
                if(pp[ad_i > 0 & ad_i <= ad_min_for_impurity,.N] > max_impurity_samples) next
                
            }
            # "Candidate mutations must be present on both the forward and reverse strand in the offspring (ADF, ADR > 0)."
            if(!(cdt[[1,adrcol]] > 0 && cdt[[1,adfcol]] > 0)) next
            # focal child GQ and DP filter
            if(!((cdt[1,igq] >= MIN_GQ) && (cdt[1, idp] %between% c(MIN_DP/2, MAX_DP/2)))) next
            # focal child ALT count filter
            if(!(cdt[[1,adcol]] >= MIN_ALT)) next # We don't half minimum alt even if in males
            # focal child pure sample
            alti <- cdt[[1, adcol]]
            ab <- alti/(alti + cdt[1,ad_ref])
            if(ab!=1) next
            mvdt_y <- rbind(mvdt_y, cdt[1,1:5])
        }
    }
}



mv_all <- rbind(mvdt_auto, mvdt_male) %>% rbind(.,mvdt_female) %>% rbind(., mvdt_y)
mv_all <- mv_all[!is.na(CHROM)]
##########################################
### cross family check
#VARIANTS_FAM_DIR="/data/home/ywang120/myData/bbsrc/idata/denovo_mutations/dmel"

variants_file <- list.files(path = VARIANTS_FAM_DIR, pattern = paste0(SPECIES, "_[0-9]*_allsample.snp.table.gz"),
                               recursive = TRUE, full.names = TRUE)
cat("Cross family variants check:\n")
cat(paste(variants_file, collapse='\n'))
cat('\n')

vdf <- data.table()

for(vv in variants_file){
    fam <- stringr::str_split(basename(vv),'_')[[1]][2]
    tmp <- fread(vv, header=T, sep='\t')
    tmp <- tmp[,1:4]
    tmp$fam <- fam
    vdf <- rbind(vdf, tmp)
}

mv_all$cross_fam_check <- 'xxx'
for(i in 1:nrow(mv_all)){
    chr = mv_all[i, CHROM]
    pos = mv_all[i, POS]
    alt = mv_all[i, ALT]
    tmp <- vdf[CHROM==chr & POS==pos]
    tmp <- tmp[fam!=FAM]
    if(nrow(tmp)!=0 ){
        # print(i)
        # break
        ifam <- paste(tmp$fam, collapse=',')
        alts <-  paste(tmp$ALT, collapse=',') %>% stringr::str_split(.,',')
        alts <- alts[[1]]
        if(alt %in% alts){
            mv_all[i,cross_fam_check:=ifam]
        }   
    }
}

mv_all <- mv_all[cross_fam_check=='xxx', 1:5]

outdf <- left_join(mv_all, aa, by = c('CHROM', 'POS', 'REF', 'ALT'))
fwrite(outdf, outfile, col.names=T, row.names=F, sep="\t", quote=F)


