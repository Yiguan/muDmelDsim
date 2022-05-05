#!/usr/bin/env Rscript


library("data.table")
library("dplyr")


sams <- list.files(pattern='sample_.*.minor.sam')

myfread <- function(ff){
    Lines <- fread(ff, sep = "", header=F)[[1]]
    n <- max(count.fields(textConnection(Lines), sep = "\t", comment.char=''))
    hh <- paste0(rep('V',n), 1:n)
    fsam <- fread(text = c(paste0(hh, collapse = '\t'), Lines), header = TRUE, fill = TRUE, sep='\t')
    return(fsam)
}


for(sam in sams){
    sid <- stringr::str_split(sam, '\\.')[[1]][1]
    outname <- paste0(sid,'.minor.mutated.sam')
    print(sam)
    fsam <- myfread(sam)
    mutated_sam <- paste0('merged_', sid,'_mutated.sam')
    fmut <- myfread(mutated_sam)
    fmut <- fmut[,c(1:6,10)]
    names(fmut)[7] <- 'H10'
    fmerge <- left_join(fsam, fmut, by =c('V1','V2','V3','V4', 'V5', 'V6')) %>% as.data.table
    fmerge[!is.na(H10),V10:=H10]
    fmerge$H10 <- NULL
    fmerge <- fmerge %>% mutate_if(is.integer, as.character) ### avoid unknown space before integer output!!!!!
    
    # write.table(stringr::str_trim(apply(fmerge, 1, paste, collapse='\t')),
    #         outname,
    #         row.names=FALSE, col.names=F, quote=F)
    writeLines(stringr::str_trim(apply(fmerge, 1, paste, collapse='\t')),
        outname)
    #fwrite(fmerge, outname, col.names=F, row.names=F, quote=F, sep='\t')
}