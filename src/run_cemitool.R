run_cemitool <- function(gse_id, expr, annot, results_dir, p){
    
    p$expr <- expr
    p$annot <- annot
    p$plot_diagnostics <- FALSE
    
    possible_error <- tryCatch(cem <- do.call(cemitool, p),
                               error=function(e) e)
    if(inherits(possible_error, "error")){
        if(!dir.exists(file.path(results_dir, gse_id))){
            dir.create(file.path(results_dir, gse_id))
        }
        cat(paste(gse_id, as.character(possible_error)), 
            file=file.path(results_dir, gse_id, "log.txt"), append=TRUE)
    }else{
        cem <- possible_error
        possible_error2 <- tryCatch({
            generate_report(cem, directory=file.path(results_dir, gse_id, "Reports"))
            write_files(cem, directory=file.path(results_dir, gse_id, "Tables"))
            save_plots(cem, directory=file.path(results_dir, gse_id, "Plots"))
            if(length(dev.list()) > 0) dev.off()
        },
        error=function(e2) e2)
        if(inherits(possible_error2, "error")){
            if(!dir.exists(file.path(results_dir, gse_id))){
                dir.create(file.path(results_dir, gse_id))
            } 
            cat(paste(gse_id, as.character(possible_error2)), 
                file=file.path(results_dir, gse_id, "log.txt"), append=TRUE)
        }else{
            return(cem)
        }
    }
}
