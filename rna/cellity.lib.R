
my_extract_features <- function (counts_nm, read_metrics, prefix = "", output_dir = "",
    common_features = NULL, GO_terms = NULL, extra_genes = NULL,
    organism = "mouse")
{
    feature_info <- get("feature_info")
    if (is.null(common_features)) {
        common_features <- feature_info[[2]]
    }
    if (is.null(GO_terms)) {
        GO_terms <- feature_info[[1]]
    }
    if (organism == "human" || organism == "org.Hs.eg.db") {
        organism <- "org.Hs.eg.db"
    }
    else if (organism == "mouse" || organism == "org.Mm.eg.db") {
        organism <- "org.Mm.eg.db"
    }
    else {
        print.warnings("You have specified a different organism than mouse or human.\n\n           This might work, but you have to make sure you have specified the appropiate database as or
ganism (e.g. org.Hs.eg.db), and you have it also installed.\n \n           Also, pleae note that extra_genes need to match the organism of interest.")
    }
    if (is.null(extra_genes)) {
        if (organism == "org.Hs.eg.db") {
            extra_genes <- get("extra_human_genes")
        }
        else if (organism == "org.Mm.eg.db") {
            extra_genes <- get("extra_mouse_genes")
        }
	}
    #.info("Extracting features")
    genes <- rownames(counts_nm)
    if (is.null(genes) | length(genes) == 0) {
        #.info("Please annotate your expression matrix with genes identifiers as rownames")
        cat("Please annotate your expression matrix with genes identifiers as rownames\n")
        return(NULL)
    }
    features_all <<- cellity:::feature_generation(counts_nm, read_metrics,
        GO_terms, extra_genes, organism)
    #.info(paste0("Features extracted."))
    cat(paste0("Features extracted.\n"))
    sds <<- apply(features_all, 2, sd)
    features_all <- features_all[, which(sds != 0)]
    types <- c("all", "common")
    features_common <- features_all[, which(colnames(features_all) %in%
        common_features)]
    if (prefix != "" && output_dir != "") {
        o <- paste(output_dir, prefix, sep = "/")
        print(output_dir)
        dir.create(o, showWarnings = FALSE, recursive = TRUE)
        f_all <- file.path(o, paste0(prefix, ".", types[1], ".features"))
        f_common <- file.path(o, paste0(prefix, ".", types[2],
            ".features"))
        write.table(features_common, f_common,sep = "\t",quote = F)
        write.table(features_all, f_all,sep = "\t",quote = F)
        cat(paste0("Features saved: \n", f_all))
        cat(paste0("Features saved: \n", f_common))
    }
    return(list(features_all, features_common))
}
		

  
