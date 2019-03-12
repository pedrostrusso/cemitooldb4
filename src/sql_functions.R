get_user_data <- function(user, dbname, host){
    password <- getPass::getPass()
    user_data <- list(user, password, dbname, host)
    names(user_data) <- c("user", "password", "dbname", "host")
    return(user_data)
}

truncate_tables <- function(pool){
    rs <- pool::dbExecute(pool, "SET foreign_key_checks = 0;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE study_info;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE cemitool_run;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE keywords;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE sample_annot;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE modules;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE ora;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE module_genes;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE enrichment;")
    rs <- pool::dbExecute(pool, "TRUNCATE TABLE interactions;")
    rs <- pool::dbExecute(pool, "SET foreign_key_checks = 1;")
}


