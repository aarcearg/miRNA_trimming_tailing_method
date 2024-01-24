f_load_trim_tail_RData <- function(
   filepath ## file path with all RData files to analyze
   ,dir_pattern="trim_tail.RData" ## pattern to search for RData files in the directory
   ,study_name_rm="_trim_tail.RData" ## text in the RData filename to remove, the rest will be used as Study name
   )
{
    ## Load trim tail tables in RData files (for one or more studies in a directory). These files correspond to the parameter output_fn in count_trim_tail.R
    library(dplyr)
    fnv <- dir(filepath,pattern=dir_pattern)
    studyv <- gsub(study_name_rm,"",fnv)
    (load(paste0(filepath,studyv[1],study_name_rm))) ## loading first RData file
    colnames(dfr_trim_tail)[4] <- studyv[1]
    dfr_trim_tail <- dfr_trim_tail[,1:4]
    dfr_trim_tail_totals <- dfr_trim_tail
    for(study in studyv[-1]){
        (load(paste0(filepath,study,study_name_rm))) ## loading RData file
        dfr_trim_tail <- dfr_trim_tail[,1:4]
        colnames(dfr_trim_tail)[4] <- study
        dfr_trim_tail_totals <- dfr_trim_tail_totals %>%
            left_join(dfr_trim_tail)
    }
    dfr_trim_tail_totals <- dfr_trim_tail_totals
    dfr_trim_tail_totals$Trimming <- -dfr_trim_tail_totals$Trimming
    return(dfr_trim_tail_totals) ## single table with all results
}


f_plot_trim_tail <- function(
    filepath="~/Dropbox/" ## directory where to store the plot
   ,filename="trimtail_globals_proportion" ## plot filename without the pdf extension (which is added)
   ,plotwidth=12,plotheight=12
   ,dfr_trim_tail_totals ## table from f_load_trim_tail_RData
   ,add_guide_star_info = FALSE ## add guide star info
   ,miRoksinStar.collapsed=NULL ## vector of GUIDE miRNAs. Provide if adding guide/star label
   ,replace_replicate=TRUE
   ,linerep_split="r"
   ,lines2plot # vector of samples to be plot
   ,do_facet_line_replicate=FALSE
   ,do_Line_relevel=FALSE
   ,ordered_line_levels=NULL ## reorder Lines
   ,filter_total_sum=FALSE
   ,total_sum_cut=NULL
   ,max_size_area=10
    )
{
    library("gridExtra")
    pdf(paste0(filepath,filename,".pdf")
       ,plotwidth,plotheight)
    miRuv <- unique(dfr_trim_tail_totals$miRNA)
    for(i in seq_along(miRuv)){
        miR <- miRuv[i]
        isstar <- ifelse(miR %in% miRoksinStar.collapsed,"GUIDE", "STAR") 
        if(i%%50==1) cat(miR," - ",i," miRNA\n")
        dfr_trim_tail_1miR <- dfr_trim_tail_totals %>%
            filter(miRNA==miR)
        if(replace_replicate) colnames(dfr_trim_tail_1miR) <- gsub("-","r",colnames(dfr_trim_tail_1miR))
        p <- dfr_trim_tail_1miR %>%
            pivot_longer(-c(miRNA,Tailing,Trimming),names_to="Sample",values_to="Counts")
        if(do_facet_line_replicate){            
            p <- p %>% 
                mutate(Line=sapply(Sample,function(x)strsplit(x,linerep_split)[[1]][1])
                      ,Replicate=sapply(Sample,function(x)strsplit(x,linerep_split)[[1]][2]))
            p <- p %>% filter(Line%in%lines2plot)
        }
        pp <- p %>% group_by(Sample) %>% summarize_each(TotalSum=Counts,sum)
        if(filter_total_sum){
            if(!any(pp$TotalSum > total_sum_cut)) next
        }
        p <- left_join(p,pp,by="Sample") %>%
            mutate(Proportion=Counts/TotalSum)
        if(do_Line_relevel) p$Line <- factor(p$Line,levels=ordered_line_levels,ordered=TRUE)
        ggp <- ggplot(p) +
            geom_count(aes_string(x="Trimming",y="Tailing",size="Proportion",colour="Counts")) +
            ## scale_size_continuous(limits=c(0,1)) +
            scale_size_area(limits=c(0,1),max_size=max_size_area) +
            coord_cartesian(xlim=c(-7,0.5),ylim=c(-0.5,7)) +
            labs(title=paste0(miR," - ",isstar)) +
            theme(plot.title=element_text(hjust=0.5))
        if(do_facet_line_replicate){
            ggp <- ggp + facet_grid(Replicate ~  Line,scales="free",space="free")
        }else{
            ggp <- ggp + facet_grid(. ~  Sample,scales="free",space="free")
        }
        print(ggp)
    }
    dev.off()
}

f_trimtail_index <- function(dfr_trim_tail_totals
                            ,Sample_ColName="SRR"){
    ### Calculate indices from total counts
    ptrim <- dfr_trim_tail_totals %>%
        filter(Trimming!=0) %>% 
        group_by(miRNA) %>%
        summarize_if(is.numeric,sum) %>%
        dplyr::select(-Trimming,-Tailing) %>% 
        pivot_longer(-miRNA,names_to=!!Sample_ColName,values_to="CountsTrim")
    ptail <- dfr_trim_tail_totals %>%
        filter(Tailing!=0) %>% 
        group_by(miRNA) %>%
        summarize_if(is.numeric,sum) %>% 
        dplyr::select(-Trimming,-Tailing) %>% 
        pivot_longer(-miRNA,names_to=!!Sample_ColName,values_to="CountsTail")
    ptotal <- dfr_trim_tail_totals %>%
        group_by(miRNA) %>%
        summarize_if(is.numeric,sum) %>% 
        dplyr::select(-Trimming,-Tailing) %>% 
        pivot_longer(,-miRNA,names_to=!!Sample_ColName,values_to="TotalCounts")
    trimtail_index <- full_join(ptrim,ptail) %>% full_join(ptotal) %>%
        mutate(TrimmingIndex=CountsTrim/TotalCounts
              ,TailingIndex=CountsTail/TotalCounts)
    return(trimtail_index)
}

f_trimtail_index_diff <- function(
    trimtail_index_list ## output of f_trimtail_index() with indices
   ,Sample_ColName="SRR" ## text contained in count columns (SRR is part of run accession ids) find them in the column
   ,ref_list ## each element can be a list with the name of the reference line and a vector of treatment lines
   ,runinfo ## table with run information, that will be joined to the counts table
   ,split_char="_R" ## character to spit sample information from replicate information, modify accordingly
   ,verbose=TRUE
    ){
    ### Function to compute differences of indices between samples for multiple studies
    library(dplyr)
    trimtail_index_diff_list <- list()
    for(study in names(trimtail_index_list)){
        if(verbose) cat("## ",study,"\n")
        ref_line <- ref_list[[study]]
        trimtail_index <- trimtail_index_list[[study]]
        ## runinfo_study <- runinfo %>% filter(Study==study) %>%
        trimtail_meanindex <- trimtail_index %>%
            left_join(runinfo,by="SRR") %>%
            mutate(Line=sapply(Sample,function(x)strsplit(x,split=split_char)[[1]][1])) %>%
            group_by(miRNA,Line) %>%
            summarise(MeanTrimmingIndex=mean(TrimmingIndex,na.rm=TRUE)
                     ,MeanTailingIndex=mean(TailingIndex,na.rm=TRUE))
        if(!is.list(ref_line)){
            trim_diff <- trimtail_meanindex %>% ## dfr with trimming differences
                dplyr::select(-MeanTailingIndex) %>%
                spread(key=Line,value=MeanTrimmingIndex) %>%
                mutate(across(where(is.numeric), ~ .x - !! sym(ref_line))) %>%
                dplyr::select(- !! sym(ref_line)) %>%
                gather("Line","DiffTrimmingIndex",-miRNA)
            tail_diff <- trimtail_meanindex %>% ## dfr with tailing differences
                dplyr::select(-MeanTrimmingIndex) %>%
                spread(key=Line,value=MeanTailingIndex) %>%
                mutate(across(where(is.numeric), ~ .x - !! sym(ref_line))) %>%
                dplyr::select(- !! sym(ref_line)) %>%
                gather("Line","DiffTailingIndex",-miRNA)
            trimtail_index_diff <- full_join(trim_diff,tail_diff,by=c("miRNA","Line"))
        }else{
            trimtail_index_diff_subgrp <- list()
            for(i in seq_along(ref_line)){
                rline <- names(ref_line[i])
                tline <- ref_line[[i]]
                trimtail_meanindex_subgrp <- trimtail_meanindex %>%
                    filter(Line%in%c(rline,tline))
                trim_diff <- trimtail_meanindex_subgrp %>% ## dfr with trimming differences
                    dplyr::select(-MeanTailingIndex) %>%
                    spread(key=Line,value=MeanTrimmingIndex) %>%
                    mutate(across(where(is.numeric), ~ .x - !! sym(rline))) %>%
                    dplyr::select(- !! sym(rline)) %>%
                    gather("Line","DiffTrimmingIndex",-miRNA)
                tail_diff <- trimtail_meanindex_subgrp %>% ## dfr with tailing differences
                    dplyr::select(-MeanTrimmingIndex) %>%
                    spread(key=Line,value=MeanTailingIndex) %>%
                    mutate(across(where(is.numeric), ~ .x - !! sym(rline))) %>%
                    dplyr::select(- !! sym(rline)) %>%
                    gather("Line","DiffTailingIndex",-miRNA)
                trimtail_index_diff_subgrp[[i]] <- full_join(trim_diff,tail_diff
                                                           ,by=c("miRNA","Line"))
            }
            trimtail_index_diff <- bind_rows(trimtail_index_diff_subgrp)
        }
        trimtail_index_diff_list[[study]] <- trimtail_index_diff %>% mutate(Study=study)
    }
    return(trimtail_index_diff_list)
}

f_trimtail_index_diff_test <- function(
    trimtail_index_list ## output of f_trimtail_index() with indices
   ,Sample_ColName="SRR" ## text contained in count columns (SRR is part of run accession ids) find them in the column
   ,ref_list ## each element can be a list with the name of the reference line and a vector of treatment lines ## NOT IMPLEMENTED HERE
   ,runinfo ## table with run information, that will be joined to the counts table
   ,split_char="_R" ## character to spit sample information from replicate information
   ,test_grps=c("all","guide","star")
   ,miRNA_guide ## vector of guide miRNAs, if test groups guide and star are going to be tested
   ,total_count_threshold=NULL ## set this to a number to filter miRNAs with less than this total counts
   ){
    ### Function to test differences in  mean indices between samples for multiple studies
    library(dplyr)
    library("broom")
    trimtail_index_diff_test_list <- list()
    for(i in test_grps){
        cat("### Test grp: ",i,"\n")
        trimtail_index_diff_test_list[[i]] <- list()
        for(study in names(trimtail_index_list)){
            cat(" ## study: ",study,"\n")
            ref_line <- ref_list[[study]]
            trimtail_index <- trimtail_index_list[[study]]
            if(i=="guide"){
                trimtail_index <- trimtail_index %>%
                    filter(miRNA %in% miRNA_guide)
            }else if(i=="star"){
                trimtail_index <- trimtail_index %>%
                    filter(!miRNA %in% miRNA_guide)
            }
            if(!is.null(total_count_threshold)){
                trimtail_index_mintotcnt <- trimtail_index %>%
                    group_by(miRNA) %>%
                    summarise(MinTotalCounts=min(TotalCounts)) %>%
                    filter(MinTotalCounts>=total_count_threshold)
                trimtail_index <- trimtail_index %>%
                    filter(miRNA%in%unlist(trimtail_index_mintotcnt$miRNA))
                cat("  # miRNAs above threshold: ",nrow(trimtail_index_mintotcnt),"\n")
            }
            trimtail_index_runinfo <- trimtail_index %>%
                left_join(runinfo,by="SRR") %>%
                mutate(Line=sapply(Sample,function(x)strsplit(x,split=split_char)[[1]][1])) %>%
                group_by(miRNA,Line) %>%
                summarise(MeanTrimmingIndex=mean(TrimmingIndex,na.rm=TRUE)
                         ,MeanTailingIndex=mean(TailingIndex,na.rm=TRUE))
            if(!is.list(ref_line)){
                trimtail_index_runinfo_ref <- trimtail_index_runinfo %>%
                    filter(Line==ref_line)
                trimtail_index_runinfo_treat <- trimtail_index_runinfo %>%
                    filter(Line!=ref_line)
                test_res <- full_join(trimtail_index_runinfo_ref,trimtail_index_runinfo_treat
                                     ,by=c("miRNA"),suffix=c("_ref","_treat")) %>%
                    nest(data=-Line_treat) %>% ## por qué -Line_treat
                    mutate(ttest_trimming=map(data, ~ t.test(.x$MeanTrimmingIndex_treat, .x$MeanTrimmingIndex_ref,paired=TRUE))
                          ,ttest_tailing=map(data, ~ t.test(.x$MeanTailingIndex_treat, .x$MeanTailingIndex_ref,paired=TRUE))
                          ,tidied_trimming=map(ttest_trimming,tidy)
                          ,tidied_tailing=map(ttest_tailing,tidy)
                          ,miRNAanalyzed=list(
                               str_c(unique(trimtail_index$miRNA),collapse="|")
                           )
                           )##  %>%
                    ## unnest(tidied)
            }else{
                ## break
                trimtail_index_diff_subgrp <- list()
                for(irl in seq_along(ref_line)){
                    rline <- names(ref_line[irl])
                    tline <- ref_line[[irl]]
                    trimtail_index_runinfo_ref <- trimtail_index_runinfo %>%
                        filter(Line==rline)
                    trimtail_index_runinfo_treat <- trimtail_index_runinfo %>%
                        filter(Line%in%tline)
                    trimtail_index_diff_subgrp[[irl]] <- full_join(trimtail_index_runinfo_ref,trimtail_index_runinfo_treat
                                         ,by=c("miRNA"),suffix=c("_ref","_treat")) %>%
                        nest(data=-Line_treat) %>% ## por qué -Line_treat
                        mutate(ttest_trimming=map(data, ~ t.test(.x$MeanTrimmingIndex_treat, .x$MeanTrimmingIndex_ref,paired=TRUE))
                              ,ttest_tailing=map(data, ~ t.test(.x$MeanTailingIndex_treat, .x$MeanTailingIndex_ref,paired=TRUE))
                              ,tidied_trimming=map(ttest_trimming,tidy)
                              ,tidied_tailing=map(ttest_tailing,tidy)
                              ,miRNAanalyzed=list(
                                   str_c(unique(trimtail_index$miRNA),collapse="|")
                               )
                               )##  %>%
                }
                test_res <- bind_rows(trimtail_index_diff_subgrp)
            }
            trimtail_index_diff_test_list[[i]][[study]] <- test_res %>%
                mutate(Study=study)
        }
    }
    return(trimtail_index_diff_test_list)
}
