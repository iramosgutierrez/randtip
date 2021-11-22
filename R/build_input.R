#' Function to create info given a species vector(column)
#' @export

build.info<- function(species, tree=NULL, find.ranks=TRUE, db="ncbi", mode="backbone",
                      interactive=F,genus=F){


  if(is.data.frame(species)){
  if(ncol(species)!=1){stop("Species must be provided as a character vector or single-column dataframe.")}else{species<- species[,1]}}
  species<- as.vector(species)
  if(!(is.vector(species) )){stop("Species must be provided as a character vector or single-column dataframe.")}
  if(is.list(species) ){stop("Species must be provided as a character vector or single-column dataframe.")}
  if(length(species[duplicated(species)])==1){stop("Taxon ", paste0(species[duplicated(species)], collapse=", "),  " is duplicated.")}
  if(length(species[duplicated(species)])>1 ){stop("Taxa ",  paste0(species[duplicated(species)], collapse=", "), " are duplicated.")}
  if(is.null(tree) & mode=="backbone"){stop("Parameter 'mode' is set to \"backbone\", but the backbone tree is missing. Please, provide a backbone tree.")}
  if (!is.null(tree) & !inherits(tree, "phylo")) { stop("Backbone tree must be an object of class \"phylo\".")}
  if(!(mode %in% c("list", "backbone"))) {stop("Parameter 'mode' must be \"list\" or \"backbone\" ")}



  names_df <- c("taxon", randtip::randtip_ranks(),
               "taxon1", "taxon2","rand.type","polyphyly.scheme", "use.paraphyletic",
               "use.singleton","use.stem","respect.mono", "respect.para","clump.puts",
               "prob","keep.tip")

  species<- gsub(" ", "_", species)
  if(any(substr(species, 1, 1)=="_")){
    species[which(substr(species, 1, 1)=="_")]<- substr(species[which(substr(species, 1, 1)=="_")], 2, nchar(species[which(substr(species, 1, 1)=="_")]))
  }
  if(any(substr(species, nchar(species), nchar(species))=="_")){
    species[which(substr(species, nchar(species), nchar(species))=="_")]<-
      substr(species[which(substr(species, nchar(species), nchar(species))=="_")],
             1, nchar(species[which(substr(species, nchar(species), nchar(species))=="_")])-1)

  }
  tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
  tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

  spp.in.tree<- tree$tip.label
  spp.original<- species

  if(mode=="backbone"){
        species<- c(species, spp.in.tree[!(spp.in.tree%in%species)])
  }

  info<- as.data.frame(matrix(nrow = length(species), ncol = length(names_df)))
  names(info)<- names_df



  info$taxon<- species
  onlygenus<-which(sapply(strsplit(species, "_"), length)==1)
  if(length(onlygenus)>0){
    for(t in onlygenus){
      if(paste0(info$taxon[t], "_sp.")%in% tree$tip.label){
        info$taxon[t]<-paste0(info$taxon[t], "_sp2.")}else{info$taxon[t]<-paste0(info$taxon[t], "_sp.")}

    }
  }


  if(isTRUE(genus)){
    if(any(sapply(strsplit(species, "_"), length)>1)){stop("Taxa must specify only genera for \"genus\" mode")}
    if(any(sapply(strsplit(tree$tip.label, "_"), length)>1)){stop("Tree tips must represent only genera for \"genus\" mode")}
    info$taxon<- randtip::firstword(info$taxon)
  }

  info$genus<- randtip::firstword(species)


  genera<- unique(randtip::firstword(species))

  searching.categories<- randtip::randtip_ranks()[-1]

    if(find.ranks){for(i in 1:length(genera)){
      tryCatch({
        if(interactive){
          search <- suppressMessages(taxize::classification(as.character(genera[i]),
                                                            db = db))[[1]]
        }else{
        out<-capture.output(suppressMessages(
          search <- taxize::classification(as.character(genera[i]),
                                                          db = db, rows=Inf)[[1]]))}

        for(cat in searching.categories){
          if(length(search[which(search$rank==cat), "name"])==0){info[info$genus==genera[i], cat]<-NA}else{
            cats<-search[which(search$rank==cat), "name"]
            if(length(cats)>1){cats<-cats[1]}
            info[info$genus==genera[i], cat]<- cats}
        }


      }, error=function(e){
        # Assign NA to fetching errors
        info[info$genus==genera[i], searching.categories] <- NA
      })

      # Avoid ip blocks. Taxize allows only 3 searches per second.
      Sys.sleep(0.33)

      if(!interactive){

        if(i==1){
          cat(paste0("Retrieving taxonomic information from ", db, " database.\n",
                     "0%       25%       50%       75%       100%", "\n",
                     "|---------|---------|---------|---------|", "\n")) }

          v<- seq(from=0, to=40, by=40/length(genera))
          v<-ceiling(v)
          v<- diff(v)
          cat(strrep("*", times=v[i]))

          if(i ==length(genera)){cat( "*\n")}


      }
    }}

  info[!(species%in%spp.original),
      c("taxon1", "taxon2","rand.type", "polyphyly.scheme","use.paraphyletic", "use.singleton",
        "use.stem","respect.mono","respect.para","clump.puts", "prob" )]<-"-"
  info$keep.tip[!(species%in%spp.original)]<- 0
  info$keep.tip[  species%in%spp.original ]<- 1
  return(info)
}

#' @export
info2input<- function(info, tree){

  input<-info
  input[is.na(input$keep.tip), "keep.tip"]<-"1"
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  input<- randtip::correct.DF(input)
  input$taxon <- gsub(" ", "_", input$taxon)

  tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
  tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)

  if(any(!(notNA(info$rand.type)%in% c("random", "polytomy", "-")))){
    errortaxa<- info[!is.na(info$rand.type),]
    errortaxa<- errortaxa$taxon[!(errortaxa$rand.type%in%c("random", "polytomy", "-"))]
    stop(paste0("\"rand.type\" argument must be \'random\' or \'polytomy\'. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa, collapse = "\n")))
  }
  if(any(!(notNA(info$polyphyly.scheme)%in% c("complete", "largest","frequentist", "-")))){
    errortaxa<- info[!is.na(info$polyphyly.scheme),]
    errortaxa<- errortaxa$taxon[!(errortaxa$polyphyly.scheme%in%c("complete", "largest","frequentist", "-"))]
    stop(paste0("\"polyphyly.scheme\" argument must be \'largest\', \'frequentist\' or \'complete\'. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa, collapse =  "\n")))
  }
  if(any(!(notNA(info$use.paraphyletic)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$use.paraphyletic),]
    errortaxa<- errortaxa$taxon[!(errortaxa$use.paraphyletic%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"use.paraphyletic\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa,  collapse =  "\n")))
  }
  if(any(!(notNA(info$use.singleton)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$use.singleton),]
    errortaxa<- errortaxa$taxon[!(errortaxa$use.singleton%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"use.singleton\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa, collapse = "\n")))
  }
  if(any(!(notNA(info$use.stem)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$use.stem),]
    errortaxa<- errortaxa$taxon[!(errortaxa$use.stem%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"use.stem\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa,  collapse = "\n")))
  }
  if(any(!(notNA(info$respect.mono)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$respect.mono),]
    errortaxa<- errortaxa$taxon[!(errortaxa$respect.mono%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"respect.mono\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa,  collapse = "\n")))
  }
  if(any(!(notNA(info$respect.para)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$respect.para),]
    errortaxa<- errortaxa$taxon[!(errortaxa$respect.para%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"respect.para\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa,  collapse = "\n")))
  }
  if(any(!(notNA(info$clump.puts)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$clump.puts),]
    errortaxa<- errortaxa$taxon[!(errortaxa$clump.puts%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"clump.puts\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa, collapse =  "\n")))
  }
  if(any(!(notNA(info$prob)%in% c("TRUE", "FALSE", "-")))){
    errortaxa<- info[!is.na(info$prob),]
    errortaxa<- errortaxa$taxon[!(errortaxa$prob%in%c("TRUE", "FALSE", "-"))]
    stop(paste0("\"prob\" argument must be logical. ",
                "Please check your info at the following taxa:\n",
                paste0(errortaxa,  collapse = "\n")))
  }

  input$MDCC     <- as.character(NA)
  input$MDCC.rank <- as.character(NA)

  input$MDCC[input$taxon%in%tree$tip.label]     <- "Tip"
  input$MDCC.rank[input$taxon%in%tree$tip.label] <- "Tip"


  input_search<- randtip::usingMDCCfinder(input = input, taxon = input$taxon[!(input$taxon%in%tree$tip.label)], tree = tree)

  input$MDCC[!(input$taxon%in%tree$tip.label)]     <- input_search[[1]]
  input$MDCC.rank[!(input$taxon%in%tree$tip.label)] <- input_search[[2]]

  #3.2 Taxa with no MDCC

  not.included<- input[is.na(input$MDCC),]
  if(length(not.included$taxon) > 0){
    message("The following taxa were not assigned MDCC and will not be bound to the tree:\n",
            paste0(not.included$taxon, "\n"))}

  return(input)
}



