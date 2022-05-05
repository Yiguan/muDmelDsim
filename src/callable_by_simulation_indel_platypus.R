#!/usr/bin/env Rscript
library("data.table")
library("dplyr")
library("tidyr")
library("argparse")

parser <- ArgumentParser()
parser$add_argument("--infile", help="input file with relevant annotation from GATK, *table.gz")
parser$add_argument("--pedfile", help="pedigree file, all samples from one family!")
parser$add_argument("--chrfile", help="comma separated file, telling whether auto or X chromosome")
parser$add_argument("--ref_fai", help="reference genome index file")
parser$add_argument("--mutation_table", help="The original table of all random generated sites")
parser$add_argument("--mutation_results", help="A table whether mutations were successfully created")
parser$add_argument("--recover_file", help='Result file after running denovo filtering on simulated data.')
parser$add_argument("--simSummary_outfile", help="A table summary the outcomes for all simulated mutations")
parser$add_argument("--callable_outfile", help="ouput file with condidate donovo")

args <- parser$parse_args()

pedfile = args$pedfile #'dmel_fam1.ped'
infile = args$infile #/sim_dmel_27_platypus_indels_biallelic.table.gz'
fai_file = args$ref_fai #'mel-all-chromosome-r5.44.fasta.fai'
sim_file = args$mutation_table #'dmel_27_random_site_mutations_indel.txt'
sim_result = args$mutation_results # 'dmel_27_modifySAM_snp.results'
recover_file = args$recover_file # 'sim_dmel_27_platypus_indels_biallelic_denovo.txt'
chrfile = args$chrfile
simSummary_outfile = args$simSummary_outfile
callable_outfile = args$callable_outfile
##########################################################
MIN_DP = 10
MAX_DP = 150
START_EXTEND = 5 

# SPECIES = stringr::str_split(basename(pedfile),'_')[[1]][1]
# FAM = stringr::str_split(stringr::str_split(basename(pedfile),'_')[[1]][2],'\\.')[[1]][1]



#############################################################
Rcpp::cppFunction("
Rcpp::IntegerVector adjustPOS(Rcpp::StringVector sid, Rcpp::StringVector chr, Rcpp::IntegerVector pos_low,
                             Rcpp::IntegerVector pos_up, Rcpp::IntegerVector pos,
                             Rcpp::StringVector SID, Rcpp::StringVector CHR, Rcpp::IntegerVector POS)
{
    int n_size = sid.size();
    int N_SIZE = SID.size();
    Rcpp::IntegerVector RET(N_SIZE, 0);
    for(int i=0; i < N_SIZE; i++)
    {
        for(int j=0; j<n_size; j++)
        {
            if(SID[i]==sid[j] && CHR[i]==chr[j] &&  POS[i] >= pos_low[j] && POS[i] <= pos_up[j])
            {
                RET[i] = pos[j];
                break;
            }
        }
    }
    return RET;
}
"
)
###############################################################
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

###################################################
aa <- fread(sim_file, header=T, sep='\t')
names(aa) <- c('INDEL_TYPE', 'INDEL_LEN','SID','SNPIDX', 'CHROM', 'POS')
###
bb <- fread(sim_result, header=T, sep='\t')
bb <- bb %>% separate('file', c('t1','t2','t3'), sep="_") %>% 
    separate('t3', c('t5','t6'), sep="\\.") %>% 
    mutate(t6=NULL) %>%
    mutate(t5 = as.numeric(t5)) %>% 
    setNames(c('SID','CHROM', 'POS','result'))
ab <- left_join(aa, bb, by = c('SID','CHROM','POS'))
##

#recover_snp <- fread(recover_file, header=T, sep = '\t')

#############################################################


cc <- fread(infile, header = TRUE, sep = "\t")
names(cc) <- stringr::str_replace_all(names(cc), ms, "mother")
names(cc) <- stringr::str_replace_all(names(cc), fs, "father")
cc <- cc[!ALT %like% ',']
called <- cc[,c('CHROM','POS')]

ab$pos_low <- ab$POS - START_EXTEND
ab$pos_up <- ab$POS + START_EXTEND

ab$called <- 0

for(i in 1:nrow(ab)){
    ichr = ab[i,CHROM]
    ipos_low = ab[i,pos_low]
    ipos_up = ab[i, pos_up]
    if(called[CHROM==ichr & POS %between% c(ipos_low, ipos_up),.N] >= 1){
        ab[i,called:=1]
    }
}

# Check parental depth and purity
ab$parentalOK <- -9

c_auto <- cc[CHROM %in% AUTO_CHR & 
            mother.NR %between% c(MIN_DP, MAX_DP)  & mother.NV==0 &
            father.NR %between% c(MIN_DP, MAX_DP)  & father.NV==0]

c_x <- cc[CHROM %in% X_CHR & mother.NV==0 & father.NV==0]
c_x_male <- c_x[mother.NR %between% c(MIN_DP, MAX_DP) ]
c_x_female <- c_x[mother.NR %between% c(MIN_DP, MAX_DP)  &
                 father.NR  %between% c(MIN_DP/2, MAX_DP/2)]


for(i in 1:nrow(ab)){
    if(ab[i,result] != 'Succeed') {next}
    if(ab[i,called]==0){next}
    ichr <- ab[i,CHROM]
    sid <- ab[i, SID]
    pos_low <- ab[i, pos_low]
    pos_up <- ab[i, pos_up]
    if(ichr %in% AUTO_CHR){
        if(c_auto[CHROM==ichr & POS %between% c(pos_low, pos_up),.N]>=1){
            ab[i, parentalOK:=1]
        }else {
            ab[i, parentalOK:=0]
        }
    }else if (ichr %in% X_CHR) {
        if(sid %in% male_offspring){
            if(c_x_male[CHROM==ichr & POS %between% c(pos_low, pos_up), .N]>=1){
                ab[i, parentalOK:=1]
            }else {
                ab[i, parentalOK:=0]
            }
        } else if (sid %in% female_offspring) {
            if(c_x_female[CHROM==ichr & POS %between% c(pos_low, pos_up),.N]>=1){
                ab[i, parentalOK:=1]
            }else {
                ab[i, parentalOK:=0]
            }
        } else {
            cat("Cannot determine sample gender: ", sid, "\n")
        }
    } else {
        cat("Cannot determine autosome or sexual chromosome: ", ichr, "\n")
    }  
}

################################################

dd <- fread(recover_file, header = T, sep = '\t')
dd$SID <- dd %>% select(ends_with('NV')) %>% apply(.,1, function(x) names(.)[which(x!=0)]) %>%
    stringr::str_replace_all('.NV$','')

dd <- dd[,c('CHROM', 'POS', 'REF', 'ALT', 'SID')]

dd[, INDEL_LEN:=abs(nchar(REF) - nchar(ALT))]
dd[nchar(REF)-nchar(ALT)>0, INDEL_TYPE:='deletion']
dd[nchar(REF)-nchar(ALT)<0, INDEL_TYPE:='insertion']

dd$adjust_pos <- adjustPOS(ab$SID, ab$CHROM, ab$pos_low, ab$pos_up, ab$POS, dd$SID, dd$CHROM, dd$POS)

dd <- dd[adjust_pos ==0 , adjust_pos:=POS]

dd <- dd[,c('SID','INDEL_TYPE','CHROM','adjust_pos')] %>% unique()
dd$detected <- 1



ab$adjust_pos <- ab$POS

abcd <- left_join(ab, dd, by=c('SID', 'INDEL_TYPE','CHROM','adjust_pos'))
discovery_rate <- abcd[result=='Succeed' & called==1 & parentalOK==1 & !is.na(detected),.N]/abcd[result=='Succeed' & called==1 & parentalOK==1,.N]
cat('\n===============\n')
cat('Recovery rate for the simulations:\n')
cat(discovery_rate)
cat('\n===============\n')
abcd[is.na(detected), detected:=0]
fwrite(abcd, simSummary_outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')

##########################
###### Callable estimate

sim_callable <- abcd[result=='Succeed' & called==1 & parentalOK==1 ]

# callable proportion
auto_call_prop <- sim_callable[CHROM %in% AUTO_CHR,.N]/aa[CHROM %in% AUTO_CHR,.N]
x_male_prop <- sim_callable[CHROM %in% X_CHR & SID %in% male_offspring,.N]/aa[CHROM %in% X_CHR & SID %in% male_offspring,.N]
x_female_prop <- sim_callable[CHROM %in% X_CHR & SID %in% female_offspring,.N]/aa[CHROM %in% X_CHR & SID %in% female_offspring,.N]
# callable count

offspring <- unique(aa$SID)
auto_numbers <- length(offspring)
auto_len <- sum(fai[V1 %in% AUTO_CHR, V2])
auto_count <- auto_numbers * 2 * auto_len * auto_call_prop
x_len <- sum(fai[V1%in% X_CHR, V2])
male_numbers <- length(offspring[offspring %in% male_offspring])
x_male_count <-  male_numbers * 1 * x_len * x_male_prop
female_numbers <- length(offspring[offspring %in% female_offspring])
x_female_count <- female_numbers * 2 * x_len * x_female_prop

total_count <- auto_count + ifelse(is.nan(x_male_count),0,x_male_count) + 
    ifelse(is.nan(x_female_count),0,x_female_count)

cat('\n===============\n')
cat('The total callable sites:\n')
cat(total_count)
cat('\n===============\n')
# out format
out_chr = c('autosome','X_male','X_female')
chr_len = c(auto_len, x_len, x_len)
sample_n = c(auto_numbers, male_numbers, female_numbers)
call_prop <- c(auto_call_prop, x_male_prop, x_female_prop)
total_call <- c(auto_count, x_male_count, x_female_count)
outdf <- data.table('chr_class' = out_chr, 'ref_len' = chr_len, 'sample_n' = sample_n,
                    'call_prop' = call_prop, 'total_called' = total_call)
fwrite(outdf, callable_outfile, col.names=TRUE, row.names=FALSE, quote=FALSE, sep = '\t')



##########################
###### Plot venn

callable <- abcd[result=='Succeed' & called==1 & parentalOK==1 ] %>% 
    tidyr::unite('uid', c('SID', 'CHROM','adjust_pos','INDEL_TYPE'))

tmp <- dd %>%  tidyr::unite('uid', c('SID','CHROM','adjust_pos', 'INDEL_TYPE'))

library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

png('simulation_venn.png')
grid.newpage()
venn.plot <- draw.pairwise.venn(
	area1 = length(tmp$uid), 
	area2 = length(callable$uid),
	cross.area = length(intersect(tmp$uid, callable$uid)),
	category = c("Detected", "Callable"),
	fill = myCol[1:2],
	cex = 1.5,
	cat.cex = 2,
	cat.pos = c(180, 180),
	ext.line.lty = "dashed",
	scaled     = FALSE,
	rotation.degree = ifelse(length(callable$uid)>length(tmp$uid),180,0)
	)
grid.draw(venn.plot)
dev.off()