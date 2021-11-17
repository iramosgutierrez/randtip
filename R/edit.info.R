

edit.info <- function (info, PUTs, remove=FALSE, column =NULL, value = NULL){

  info <- randtip::correct.DF(info)
  PUTs<- gsub(" ", "_", PUTs)

  if(length(PUTs[!(PUTs %in% info$taxon)])==1){
    stop("PUT ", PUTs[!(PUTs %in% info$taxon)], " is not included in info$taxon.")}
  if(length(PUTs[!(PUTs %in% info$taxon)]) >1){
    stop("PUTs ", paste0("\"",PUTs[!(PUTs %in% info$taxon)], "\"", collapse = ", "), " are not included in info$taxon.")}

  if(isTRUE(remove)){return(info[!(info$taxon%in%PUTs),])}

  if(is.null(column)|is.null(value)){stop("Both \'column\' and \'value\' arguments must be specified")}

  if(!(column %in% names(info))){stop("Specified \'column\' value is not a correct info column name.")}

  if(column=="rand.type" & !(value %in% c("random", "polytomy"))){stop("Argument 'rand.type' must be \"random\" or \"polytomy\"" )}
  if(column=="polyphyly.scheme" & !(value %in% c("complete", "largest", "frequentist"))){stop("Argument 'polyphyly.scheme' must be \"frequentist\", \"complete\" or \"largest\" ")}

  if(column=="use.paraphyletic" & !is.logical(value)){stop("Argument 'use.paraphyletic' must be logical." )}
  if(column=="use.singleton"    & !is.logical(value)){stop("Argument 'use.singleton' must be logical." )}
  if(column=="use.stem"         & !is.logical(value)){stop("Argument 'use.stem' must be logical." )}
  if(column=="respect.mono"     & !is.logical(value)){stop("Argument 'respect.mono' must be logical." )}
  if(column=="respect.para"     & !is.logical(value)){stop("Argument 'respect.para' must be logical." )}
  if(column=="clump.puts"       & !is.logical(value)){stop("Argument 'clump.puts' must be logical." )}
  if(column=="prob"             & !is.logical(value)){stop("Argument 'prob' must be logical." )}

  info[which(info$taxon %in% PUTs), column] <- value

  return(info)

  }
