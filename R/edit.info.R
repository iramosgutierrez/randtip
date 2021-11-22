

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

  if(isTRUE(remove.row)){return(info[!(info$taxon%in%PUTs),])}

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

  info[which(info$taxon %in% PUTs), column] <- edit

  return(info)

}


edit.tree <- function(tree,tips, edit=NULL, remove.tips=FALSE) {

  tips<- gsub(" ", "_", tips)


  if(length(tips[!(tips %in% tree$tip.label)])==1){
    stop("PUT ", tips[!(tips %in% tree$tip.label)], " is not included in the tree.")}
  if(length(tips[!(tips %in% tree$tip.label)]) >1){
    stop("tips ", paste0("\"",tips[!(tips %in% tree$tip.label)], "\"", collapse = ", "), " are not included in the tree.")}


  if(isTRUE(remove.tips)){tree <- ape::drop.tip(tree, tips); return(tree)}

  if(!(is.null(edit)) & length(tips)>1){stop("Tree tip editions must be performed individually.")}

  if(!(is.null(edit))){tree$tip.label[which(tree$tip.label==tips)]<- edit}
  return(tree)

}
