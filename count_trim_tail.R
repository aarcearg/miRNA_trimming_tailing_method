## Run: Rscript count_trim_tail.R read_counts.txt.gz output_fn_prefix miRNA_collapsed.fa

library(Biostrings)
library(dplyr)
args <- commandArgs(trailingOnly=TRUE) ## loading arguments

## test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("No arguments were supplied.\n", call.=FALSE)
}

run_read_counts <- args[1] ## input read counts filename

output_fn <- args[2] ## output filename prefix for the study (do not include .RData in this filename)

## reading miRNA fasta file
miR_collapsed <- readDNAStringSet(args[3])

## loading reads
con <- gzfile(run_read_counts,open="r")
p <- read.table(con)
close(con)

readsStrSet <- DNAStringSet(p[,2])
mcols(readsStrSet) <- p[,1]
colnames(mcols(readsStrSet)) <- "Counts"
rm(p)

dfr_trim_tail <- data.frame(matrix(nrow=0,ncol=8))
colnames(dfr_trim_tail) <- c("miRNA","Tailling","Trimming","Total","A","C","G","T")

## checking a "seed" sequence of 13 nt that should match the miRNA to analyze
readsStrSet <- readsStrSet[width(readsStrSet) >= 13]
readsStrSet_13seed <- subseq(readsStrSet,1,13)
miR_collapsed_13seed_unq <- unique(subseq(miR_collapsed,1,13))
idx <- which(readsStrSet_13seed %in% miR_collapsed_13seed_unq)
readsStrSet_miR13seed <- readsStrSet[idx]

cat("Total counts in readsStrSet: ",sum(mcols(readsStrSet)$Counts),"\n")
cat("Total counts in readsStrSet_miR13seed: ",sum(mcols(readsStrSet_miR13seed)$Counts),"\n")

## Creating an empty data frame to fill with counts
dfr_trim_tail <- data.frame(miRNA=rep(names(miR_collapsed),each=64)
                           ,Trimming=rep(0:7,each=8)
                           ,Tailing=rep(0:7,times=8)
                           ,Total=0
                           ,A=0,C=0,G=0,T=0
                            )

## list of trimmed miRNAs
miR_collapsed_trimmed_list <- lapply(0:7, function(x)subseq(miR_collapsed,1,-1-x))

for(i in 1:length(readsStrSet_miR13seed)){
    ## cat("Iteration: ",i,"\n")
    read <- readsStrSet_miR13seed[i]
    read_counts <- mcols(read)$Counts
    tail <- 0
    ## I start without tailing
    any_trimmed <- FALSE
    lastBase <- as.character(subseq(read,-1,-1))
    if(lastBase=="N") next
    for(trim in 0:7){
        miR_collapsed_trimmed <- miR_collapsed_trimmed_list[[trim+1]]
        idx_miR <- miR_collapsed_trimmed %in% as.character(read)
        if(any(idx_miR)){
            miRs <- names(miR_collapsed)[idx_miR]
            irow_dfr <- which(dfr_trim_tail$miRNA %in% miRs &
                              dfr_trim_tail$Trimming == trim &
                              dfr_trim_tail$Tailing == tail)
            dfr_trim_tail[irow_dfr,"Total"] <- dfr_trim_tail[irow_dfr,"Total"] + read_counts
            dfr_trim_tail[irow_dfr,lastBase] <- dfr_trim_tail[irow_dfr,lastBase] + read_counts
            any_trimmed <- TRUE
        }
        ## if(any_trimmed) break ##
    }
    if(any_trimmed) next ## If read matched with a trimmed or untrimmed miR, the read is not further analyzed
    any_hit <- FALSE
    for(trim in 0:7){
        ## cat("Trim: ",trim,"\n")
        read_trimmed_list <- lapply(1:7, function(x)subseq(read,1,-1-x))
        miR_collapsed_trimmed <- miR_collapsed_trimmed_list[[trim+1]]
        for(tail in 1:7){
            ## cat("Tail: ",tail,"\n")
            read_trimmed <- read_trimmed_list[[tail]]
            idx_miR <- miR_collapsed_trimmed %in% as.character(read_trimmed)
            if(any(idx_miR)){
                miRs <- names(miR_collapsed)[idx_miR]
                irow_dfr <- which(dfr_trim_tail$miRNA %in% miRs &
                                  dfr_trim_tail$Trimming == trim &
                                  dfr_trim_tail$Tailing == tail)
                dfr_trim_tail[irow_dfr,"Total"] <- dfr_trim_tail[irow_dfr,"Total"] + read_counts
                dfr_trim_tail[irow_dfr,lastBase] <- dfr_trim_tail[irow_dfr,lastBase] + read_counts
                any_hit <- TRUE
                break
            }
        }
        if(any_hit) break
    }
}
cat("Total counts (dfr_trim_tail$Total): ",sum(dfr_trim_tail$Total),"\n")

## the output table is saved as RData file (with the output filename provided)
save(dfr_trim_tail,file=paste0(output_fn,"trim_tail.RData")) 
