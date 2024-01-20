
#' Function to edit info data frames.
#'
#' Auxiliar function for smoothly editing the info data frame
#'
#' @param info An info data frame.
#' @param taxa Character vector with the taxa to be edited or removed.
#' @param column The name of the column to be edited for the specified taxa.
#'               NULL is only accepted if \code{remove.rows} is TRUE.
#' @param edit Any allowed value for the column to be edited. NULL is only
#'             accepted if \code{remove.rows} is TRUE.
#' @param remove.rows If TRUE, the specified taxa will be removed from the
#'                    info data frame (default is FALSE).
#'
#' @return An edited info data frame
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
#'  catspecies <- c("Lynx_lynx", "Panthera_uncia",
#' "Panthera_onca", "Felis_catus", "Puma_concolor",
#' "Lynx_canadensis", "Panthera_tigris", "Panthera_leo",
#' "Felis_silvestris")
#'
#' cats.info <- build_info(species=catspecies, tree= cats,
#'      find.ranks=TRUE, db="ncbi", mode="backbone")
#'
#' cats.info <- edit_info(cats.info, taxa= "Puma_concolor",
#' column = "subfamily", edit = "Felinae")
#'
#' @export
edit_info <- function (info, taxa, column =NULL, edit = NULL, remove.rows=FALSE){

    #if(file.exists(info)){                                             # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 
    #  if(grep(getwd(), info)==1){filedir <-  info}else{                # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 
    #    filedir <- paste0(getwd(), "/", info)                          # WE SHOULD EITHER REMOVE OR KEEP THESE LINES     
    #  }                                                                # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 

    #  cat(paste0("Reading info file from\n", filedir))                 # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 
    #  info <- read.table(info)                                         # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 
    #}                                                                  # WE SHOULD EITHER REMOVE OR KEEP THESE LINES 

    info <- correct_DF(info)

    taxa <- stringr::str_trim(taxa)
    taxa<- gsub(" ", "_", taxa)

    if(length(taxa[!(taxa %in% info$taxon)])>=1){
        stop("taxa ",
            paste0("\"",taxa[!(taxa %in% info$taxon)], "\"", collapse = ", "),
            " are not included in the info data frame")
    }
    if(remove.rows){return(info[!(info$taxon%in%taxa),])}
    if(is.null(column)|is.null(edit)){
        stop("Both \'column\' and \'edit\' arguments must be specified")
    }
    if(!(column %in% names(info))){
        stop("Specified \'column\' value is an invalid column name.")
    }

    rand.types <- c("random", "polytomy", NA)
    if(column=="rand.type" & !(edit %in% rand.types)){
        stop("Argument 'rand.type' must be \"random\" or \"polytomy\" or NA." )
    }
    polyphyly.schemes <- c("complete", "largest", "frequentist", NA)
    if(column=="polyphyly.scheme" & !(edit %in% polyphyly.schemes)){
        stop("Argument 'polyphyly.scheme' must be \"frequentist\", ",
            "\"complete\" or \"largest\" or NA.")
    }

    logical.cols <- c("use.paraphyletic", "use.singleton", "use.stem",
                      "respect.mono", "respect.para", "clump.puts", "prob")
    if(column %in% logical.cols & !is.logical(edit)){
        stop("Argument '", col, "' must be logical." )
    }

    if(column %in% c("taxon", "taxon1", "taxon2")){edit<- gsub(" ", "_", edit)}

    info[which(info$taxon %in% taxa), column] <- edit

    return(info)

}


#' Function to edit phylogenetic tip labels
#'
#' Auxiliar function for smoothly editing the tip labels of the backbone tree
#'
#' @param tree A backbone tree.
#' @param tips Character vector with the phylogenetic tips to be edited or removed.
#' @param edit A vector of the same length as \code{tips} with the new labels for
#'             the specified tips. NULL is only accepted if \code{remove.tips} is
#'             TRUE.
#' @param remove.tips If TRUE, the specified tips will be pruned
#'                    from the backbone tree (default is FALSE).
#'
#' @return An edited backbone tree.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examplesIf interactive()
#' cats <- edit_tree(cats, tips="Felis_silvestris",
#' edit= "Felis_silvestris_ssp._silvestris")
#'
#' @export
edit_tree <- function(tree, tips, edit=NULL, remove.tips=FALSE) {

   #if(file.exists(tree)){                                        # SAME AS COMMENTED ABOVE
   #  if(grep(getwd(), tree)==1){filedir <-  tree}else{           # SAME AS COMMENTED ABOVE
   #    filedir <- paste0(getwd(), "/", tree)                     # SAME AS COMMENTED ABOVE
   #  }                                                           # SAME AS COMMENTED ABOVE
   #                                                              # SAME AS COMMENTED ABOVE  
   #  cat(paste0("Reading tree file from\n", filedir))            # SAME AS COMMENTED ABOVE  
   #  tree <- ape::read.tree(tree)                                # SAME AS COMMENTED ABOVE  
   #}                                                             # SAME AS COMMENTED ABOVE  

    tips<- gsub(" ", "_", tips)
    edit<- gsub(" ", "_", edit)

    tips.in.tree <- (tips %in% tree$tip.label)
    if(length(tips[!(tips.in.tree)]) >= 1){
        stop("Tips ",
             paste0("\"",tips[!(tips.in.tree)], "\"", collapse = ", "),
             " are not included in the tree.")
    }

    if(isTRUE(remove.tips)){
        tree <- ape::drop.tip(tree, tips)
        return(tree)
    }

    if(is.null(edit)){
        stop("'edit' argument can only be NULL, if 'remove.tips' is TRUE.")
    }
    if(any(duplicated(edit))){
        stop("Duplicates detected in the list of tips to be edited.")
    }
    if(length(edit)!=length(tips)){
        stop("Arguments \'edit\' and \'tips\' must have the same length.")
    }

    for(tip in tips){
        edit.i<- edit[which(tips==tip)]
        tree$tip.label[which(tree$tip.label==tip)]<- edit.i
    }

    return(tree)
}
