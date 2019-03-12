library(CEMiTool)
library(data.table)
library(RMySQL)
library(pool)
library(plyr)
library(tidyverse)
library(rentrez)
library(annotator)

user_data <- get_user_data("prusso", "cemitooldb4", "127.0.0.1")
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- CEMiTool::read_gmt(gmt_fname)
## Get example interactions file
int_df <- read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool"))

p <- list()
p$interactions <- int_df
p$gmt <- gmt_in

# Get all platforms (and their studies' GSE accession) that are "in situ 
# oligonucleotide" or "spotted oligonucleotide" in "Homo sapiens"and that have 
# over 100 studies
gpl_ids0 <- get_gpl_ids(query_term = paste("AND (in situ oligonucleotide[PTYP]",
                                           "OR spotted oligonucleotide[PTYP])",
                                           "AND Homo sapiens[ORGN]"), 
                        min_study_num = 100)

gpl_ids <- gpl_ids0

studies_to_db <- function(user_data, 
                          gpl_ids, min_sample_num=18, max_sample_num=300, 
                          results_dir="/home/prusso/Documents/CSBL/cemitooldb4/results", 
                          verbose=TRUE, force=FALSE, 
                          num_series=1, ...){
    p <- list(...)
    p$verbose <- verbose
    
    if(verbose) message("Creating database link...")
    pool <- pool::dbPool(drv=RMySQL::MySQL(),
                         dbname=user_data$dbname,
                         host=user_data$host,
                         user=user_data$user,
                         password=user_data$password)
    #conn <- poolCheckout(pool)
    
    for(gpl_id in names(gpl_ids)){
        message(gpl_id)
        gpl_res_dir <- file.path(results_dir, gpl_id)
        for(gse_id in gpl_ids[[gpl_id]]){
            message(gse_id)
            
            if(gse_id %in% dir(gpl_res_dir)) next
            
            gse_res_dir <- file.path(gpl_res_dir, gse_id)
            create_res_dir(gse_res_dir)

            gse <- tryCatch(gse <- get_gse(gse_id),
                                       error=function(e) e)
            if(inherits(gse, "error")) {
                message("Error in GEOquery")
                next
            }
            
            expr <- get_expr_from_gse(gse, num_series=num_series)
            if(is.null(expr)){message("Error getting expr"); next}
            
            sample_num <- ncol(expr) - 1
            if(sample_num > max_sample_num | sample_num < min_sample_num){
                message(paste0("Sample size (", sample_num, ") for study ", gse_id,
                              " is out of the range defined by max_sample_num",
                              " and min_sample_num"))
                next
            }
            
            expr <- process_and_collapse(expr)
            annot_error <- tryCatch(annot <- annotate_gse(gse_id), error = function(e) e)
            if(inherits(annot_error, "error")) next
            annot <- annot %>%    
                rename(Cluster=Class, Class=label1, 
                       SampleName=Sample_geo_accession) %>%
                select(SampleName, Class)
            
            run_date <- lubridate::now()
            if(verbose) message("Running CEMiTool")
            cem <- run_cemitool(gse_id, expr, annot, gse_res_dir, p)
            cem_obj_dir <- normalizePath(file.path(gse_res_dir, "cem.rds"))
            saveRDS(cem, cem_obj_dir)
            
            if(!is.null(cem)){
                if(verbose) message("Populating MySQL database")
                populate_sql(pool=pool, gse_id=gse_id, user_data=user_data, 
                             run_date=run_date, cem=cem, cem_obj_dir=cem_obj_dir)
                message("Populated!")
            }
            
        }
        #poolReturn(conn)
    }
    
    
}




