
#' Edit an 'info' file using a simple function.
#'
#' Function to edit randtip's 'info' file.
#'
#' @param info An 'info' object.
#' @param PUTs Character vector with the PUT (or PUTs) to be edited or removed.
#' @param column Which column is to be edited for the specified PUTs. NULL value is only accepted if \code{remove.rows} is set to TRUE.
#' @param edit Any allowed value for the column of info that is to be edited. NULL value is only accepted if \code{remove.rows} is set to TRUE.
#' @param remove.rows If TRUE, the specified PUTs will be eliminated from the 'info' object. Default is FALSE.
#'
#' @return An 'info' object including the editions or deletions asked.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @export
edit.info <- function (info, PUTs, column =NULL, edit = NULL, remove.rows=FALSE){

  info <- randtip::correct.DF(info)
  PUTs<- gsub(" ", "_", PUTs)

  if(any(substr(PUTs, 1, 1)=="_")){
    PUTs[which(substr(PUTs, 1, 1)=="_")]<- substr(PUTs[which(substr(PUTs, 1, 1)=="_")], 2, nchar(PUTs[which(substr(PUTs, 1, 1)=="_")]))
  }
  if(any(substr(PUTs, nchar(PUTs), nchar(PUTs))=="_")){
    PUTs[which(substr(PUTs, nchar(PUTs), nchar(PUTs))=="_")]<-
      substr(PUTs[which(substr(PUTs, nchar(PUTs), nchar(PUTs))=="_")],
             1, nchar(PUTs[which(substr(PUTs, nchar(PUTs), nchar(PUTs))=="_")])-1)

  }




  if(length(PUTs[!(PUTs %in% info$taxon)])==1){
    stop("PUT ", PUTs[!(PUTs %in% info$taxon)], " is not included in info$taxon.")}
  if(length(PUTs[!(PUTs %in% info$taxon)]) >1){
    stop("PUTs ", paste0("\"",PUTs[!(PUTs %in% info$taxon)], "\"", collapse = ", "), " are not included in info$taxon.")}

  if(isTRUE(remove.rows)){return(info[!(info$taxon%in%PUTs),])}

  if(is.null(column)|is.null(edit)){stop("Both \'column\' and \'edit\' arguments must be specified")}

  if(!(column %in% names(info))){stop("Specified \'column\' value is not a correct info column name.")}

  if(column=="rand.type" & !(edit %in% c("random", "polytomy"))){stop("Argument 'rand.type' must be \"random\" or \"polytomy\"" )}
  if(column=="polyphyly.scheme" & !(edit %in% c("complete", "largest", "frequentist"))){stop("Argument 'polyphyly.scheme' must be \"frequentist\", \"complete\" or \"largest\" ")}

  if(column=="use.paraphyletic" & !is.logical(edit)){stop("Argument 'use.paraphyletic' must be logical." )}
  if(column=="use.singleton"    & !is.logical(edit)){stop("Argument 'use.singleton' must be logical." )}
  if(column=="use.stem"         & !is.logical(edit)){stop("Argument 'use.stem' must be logical." )}
  if(column=="respect.mono"     & !is.logical(edit)){stop("Argument 'respect.mono' must be logical." )}
  if(column=="respect.para"     & !is.logical(edit)){stop("Argument 'respect.para' must be logical." )}
  if(column=="clump.puts"       & !is.logical(edit)){stop("Argument 'clump.puts' must be logical." )}
  if(column=="prob"             & !is.logical(edit)){stop("Argument 'prob' must be logical." )}

  if(column %in% c("taxon", "taxon1", "taxon2")){edit<- gsub(" ", "_", edit)}

  info[which(info$taxon %in% PUTs), column] <- edit

  return(info)

}


#' Edit backbone tree using a simple function.
#'
#' Function to edit a backbone tree.
#'
#' @param tree A backbone tree.
#' @param tips Character vector with the phylogenetic tips to be edited or removed.
#' @param edit A vector of the same length as \code{tips}, including the new labels for the specified phylogenetic tips. NULL value is only accepted if \code{remove.tips} is set to TRUE.
#' @param remove.tips If TRUE, the specified tips will be pruned from the backbone tree. Default is FALSE.
#'
#' @return A backbone tree including the tip editions or deletions asked.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @export
edit.tree <- function(tree,tips, edit=NULL, remove.tips=FALSE) {

  tips<- gsub(" ", "_", tips)
  edit<- gsub(" ", "_", edit)

  if(length(tips[!(tips %in% tree$tip.label)])==1){
    stop("PUT ", tips[!(tips %in% tree$tip.label)], " is not included in the tree.")}
  if(length(tips[!(tips %in% tree$tip.label)]) >1){
    stop("tips ", paste0("\"",tips[!(tips %in% tree$tip.label)], "\"", collapse = ", "), " are not included in the tree.")}


  if(isTRUE(remove.tips)){tree <- ape::drop.tip(tree, tips); return(tree)}

  if(!(is.null(edit)) & any(duplicated(edit))){stop("Edited tree tips must not be identical.")}
  if(!(is.null(edit)) & length(edit)!=length(tips)){stop("Arguments \'edit\' and \'tips\' must have the same length")}
  if(!(is.null(edit))){

    for(tip in tips){
      edit.i<- edit[which(tips==tip)]
      tree$tip.label[which(tree$tip.label==tip)]<- edit.i
    }


    }
  return(tree)

}
