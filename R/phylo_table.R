#' function phylo.table returns a table which the main function can work 
#' with (editable by hand!).
#'
#' @export
#' @examples
#' #df <- phylo.table(species = checklist$match.name, Tree,  
#' #                  genera.phyleticity = tipos_generos, taxonomy = taxonomy)
phylo.table <- function(species, tree = NULL, taxonomy = NULL, 
                        genera.phyleticity = NULL,
                        sep = "_"){

    if(is.null(genera.phyleticity)){
        if(is.null(tree)){
            stop("tree object required")
        }else if(!inherits(tree, "phylo")){   
            stop("tree should be an object of class \"phylo\".")
        }
    }
    
    names_df <- c("taxon", "genus", "genus.type", "family", 
                  "order", "class", "aggregate.subspecies",
                  "relative.species","synonim.genus","sibling.genus", 
                  "forbidden.branches")
    phylo.df<- data.frame(matrix(nrow = length(species), ncol = length(names_df), 
                                 dimnames = list(NULL, names_df)))

    if(any(duplicated(species))){ 
      stop("There are duplicated species")
    }

    phylo.df$taxon <- species                                     
    phylo.df$genus <- stringr::word(species, 1, sep = sep)   
    taxa.genera <- unique(phylo.df$genus)

    if(is.null(taxonomy)){  
      taxonomy <- phylo.taxonomy(species = species) 
    } 
    tax.cols <- c("family", "order", "class")
    genus.taxonomy <- merge(phylo.df[, c("taxon", "genus")], 
                          taxonomy, id = "genus", 
                          all.x = TRUE, sort = FALSE)
    row.order <- match(phylo.df$taxon, genus.taxonomy$taxon)
    phylo.df[, tax.cols] <- genus.taxonomy[row.order, tax.cols]

    if(is.null(genera.phyleticity)){ 
    genera.phyleticity<- data.frame("genus" = taxa.genera, "type" = NA)
    for(i in 1:nrow(genera.phyleticity)){
      genera.phyleticity$type[i]<- phyleticity(tree,  
                                               genera.phyleticity$genus[i])
    }
    }
    genus.type <- merge(phylo.df, genera.phyleticity, id = "genus", 
                      all.x = TRUE, sort = FALSE)
    row.order <- match(phylo.df$taxon, genus.type$taxon)
    phylo.df$genus.type <- genus.type[row.order, "type"]

    return(phylo.df)
}

