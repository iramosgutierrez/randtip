
#' Edit an 'info' file using a simple function.
#'
#' Function to edit randtip's 'info' file.
#'
#' @param info An 'info' object.
#' @param taxa Character vector with the taxa to be edited or removed.
#' @param column The column name to be edited for the specified taxa.
#'               NULL value is only accepted if \code{remove.rows} is set to TRUE.
#' @param edit Any allowed value for the column of info that is to be edited.
#'             NULL value is only accepted if \code{remove.rows} is set to TRUE.
#' @param remove.rows If TRUE, the specified taxa will be eliminated from the
#'                    'info' object. Default is FALSE.
#'
#' @return An 'info' object including the editions or deletions asked.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @examples 
#' cats.info <- edit.info(cats.info, taxa= "Puma_concolor", column = "subfamily", edit = "Felinae")
#' 
#' @export edit.info
#' @export
edit.info <- function (info, taxa, column =NULL, edit = NULL, remove.rows=FALSE){

    info <- correct.DF(info)

    taxa <- stringr::str_trim(taxa)
    taxa<- gsub(" ", "_", taxa)

    if(length(taxa[!(taxa %in% info$taxon)])>=1){
        stop("taxa ",
            paste0("\"",taxa[!(taxa %in% info$taxon)], "\"", collapse = ", "),
            " are not included in the column taxon of info dataframe")
    }
    if(remove.rows){return(info[!(info$taxon%in%taxa),])}
    if(is.null(column)|is.null(edit)){
        stop("Both \'column\' and \'edit\' arguments must be specified")
    }
    if(!(column %in% names(info))){
        stop("Specified \'column\' value is not a correct info column name.")
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


#' Edit backbone tree using a simple function.
#'
#' Function to edit a backbone tree.
#'
#' @param tree A backbone tree.
#' @param tips Character vector with the phylogenetic tips to be edited or removed.
#' @param edit A vector of the same length as \code{tips}, including the
#'             new labels for the specified phylogenetic tips. NULL value is
#'             only accepted if \code{remove.tips} is set to TRUE.
#' @param remove.tips If TRUE, the specified tips will be pruned
#'                    from the backbone tree. Default is FALSE.
#'
#' @return A backbone tree including the requested tip editions or deletions.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#' 
#' @examples 
#' cats <- edit.tree(cats, tips="Felis_silvestris", edit= "Felis_silvestris_ssp._silvestris")
#'
#' @export edit.tree
#' @export
edit.tree <- function(tree, tips, edit=NULL, remove.tips=FALSE) {

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
        stop("Edited tree tips cannot be duplicated.")
    }
    if(length(edit)!=length(tips)){
        stop("Arguments \'edit\' and \'tips\' must have the same length")
    }

    for(tip in tips){
        edit.i<- edit[which(tips==tip)]
        tree$tip.label[which(tree$tip.label==tip)]<- edit.i
    }

    return(tree)
}
