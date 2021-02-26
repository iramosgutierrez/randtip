#' function phylo.taxonomy returns a table with family, order and 
#' class for every 
#' genus given in the species vector.
#' 
#' @export
#' @examples
#' phylo.taxonomy(species=c("Abies_pinsapo", "Pinus_pinaster", 
#'                          "Achyllea_milleifolium"))
phylo.taxonomy <- function(species, db="ncbi"){

    taxa.genera <- stringr::word(species, 1, sep = "_")
    taxa.genera <- unique(taxa.genera)

    taxonomy.df <- data.frame("genus" = taxa.genera, "family" = NA, 
                              "order" = NA, "class" = NA)
    for(i in 1:length(taxa.genera)){    
        tryCatch({
            # Search taxonomic categories
            message(paste0("Taxonomic information is being looked up in ", db ,
                           " Taxonomy Database. This may take a while!"))
            search <- taxize::classification(as.character(taxa.genera[i]), 
                                          db = db)[[1]]  
            fam <- search[search$rank == "family","name"]
            ord <- search[search$rank == "order" ,"name"]
            cla <- search[search$rank == "class" ,"name"]
          
            taxonomy.df[i, c("family", "order", "class")]<- c(fam, ord, cla)
         }, error=function(e){
            # Assign NA to fetching errors
            taxonomy.df[i, c("family", "order", "class")] <- NA 
        })

        # Avoid ip blocks. Taxize allows only 3 searches per second. 
        Sys.sleep(0.33)         
    }

    return(taxonomy.df)
}
