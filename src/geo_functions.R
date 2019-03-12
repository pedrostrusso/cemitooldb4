# Get all GPL ids
get_gpl_ids <- function(verbose=TRUE, query_term="", min_study_num=100){
    gds_search <- rentrez::entrez_search(db="gds", 
                                         term=paste0("gpl[ETYP] ", query_term), 
                                         use_history = TRUE)
    
    study_count <- gds_search$count
    
    platforms <- list()
    
    for(seq_start in seq(0, study_count, 40)){
        if(verbose) message(seq_start)
        search_res <- rentrez::entrez_search(db="gds", 
                                             term=paste0("gpl[ETYP] ", query_term), 
                                             use_history = TRUE,
                                             retmax=40, retstart=seq_start)
        res_ids <- search_res$ids
        if(verbose) print(res_ids)
        
        for(id in res_ids){
            summary_res <- rentrez::entrez_summary(db="gds", id=id)
            plat <- summary_res$accession
            stds <- unlist(strsplit(summary_res$gse, split=";"))
            stds <- paste0("GSE", stds)
            if(!missing(min_study_num)){
                if(length(stds) > min_study_num){
                    platforms[[plat]] <- stds    
                }    
            }else{
                platforms[[plat]] <- stds
            }
        }
    }
    return(platforms)
}

# Get all GSE ids from a GPL id
get_platform_gse_ids <- function(gpl_id){
    gpl <- GEOquery::getGEO(gpl_id)
    gpl_studies <- GEOquery::Meta(gpl)$series_id
    return(gpl_studies)
}

# Get GSE object from a GSE id
get_gse <- function(gse_id){
    gse <- GEOquery::getGEO(gse_id, GSEMatrix=TRUE, getGPL=TRUE)
    return(gse)
}

# Get expression from a GSE object
get_expr_from_gse <- function(gse, num_series=1){
    if(length(gse) > num_series){
        return(NULL)
    }else{
        expr <- as.data.frame(Biobase::exprs(gse[[1]]))
        if(nrow(expr)>0){
            expr$GeneSymbol <- Biobase::fData(gse[[1]])[["Gene Symbol"]]
            return(expr)
        }else{
            return(NULL)
        }
    }
}
