
#' Create a 'info' data frame.
#'
#' Function to create an 'info' object given a list of species.
#'
#' @param species A character vector or a single-column data frame including the species of interest. Word breakers must be blanks (" ") or underscores ("_").
#' @param tree A 'phylo' object backbone tree. It can be set to NULL if \code{mode} is set to "list".
#' @param find.ranks Logical. If TRUE, taxonomic information will be retrieved to identify supra-generic MDCCs for the PUTs.
#' @param db Taxonomic data base to search into if \code{find.ranks} is set to TRUE.Accepted values are 'ncbi' (default), 'itis', 'gbif' and 'bold'.
#' @param mode If mode is set to "list", the info file will be created using only the species given in the \code{species} argument.
#'  If "backbone" mode is specified, 'info' will also include all the tips included in the backbone tree.
#' @param interactive Logical. Whether or not ambiguous species names will be resolved manually by the user or filled in automatically with 'NA' when retrieving taxonomic information.
#' @param genus Logical. Whether or not a genus-level backbone tree is to be expanded. If set to TRUE, all tips in the backbone tree and taxa in the species vector must represent genera.
#'
#' @return A randtip 'info' data frame
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
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



#' Function to check randtip's files.
#'
#' Function to account for the PUT status of the species in 'info', spelling errors, putative MDCCs and the phyletic nature of groups of PPCR species.
#'
#' @param info An 'info' object.
#' @param tree The original backbone tree.
#' @param sim Name simmilarity threshold to detect possible misspellings on tip labels. Default value is 0.8.
#'
#' @return A data frame containing possible typographic errors, taxonomic ranks extracted from 'info' and the phyletic nature of each of them.
#' #'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
#' @export

check.info<- function(info, tree, sim=0.8){

  if(is.null(info)){stop("Data frame 'info' is missing.")}
  if(is.null(tree)){stop("Backbone tree is missing.")}

  info<- randtip::correct.DF(info)
  info$keep.tip[is.na(info$keep.tip)]<-"1"
  info$taxon<- gsub(" ", "_", info$taxon)

  tree$tip.label<- gsub(" ", "_", tree$tip.label)
  info.taxa<-info$taxon
  tree.taxa<- tree$tip.label

  DF<- info[info$keep.tip=="1",c("taxon",randtip::randtip_ranks())]
  DF$PUT.status<- NA
  DF$Typo<- F
  DF$Typo.names<- NA




  #1st we search included taxa
  for(i in 1:nrow(DF)){
    tax<- DF$taxon[i]
    if(tax %in%tree.taxa){DF$PUT.status[i]<- "Tip"}else{DF$PUT.status[i]<- "PUT"}
  }

  #2nd we look for name simmilarities
  for(i in 1:nrow(DF)){
    if(DF$PUT.status[i]=="PUT"){
      tax<-DF$taxon[i]
      sim.search<-tree.taxa[stringdist::stringsim(tree.taxa,tax)>sim]
      if(length(sim.search)>0){
        DF$Typo[i]<- TRUE
        sim.search<-paste0(sim.search, collapse = " / ")
        DF$Typo.names[i]<- sim.search}
      rm(tax, sim.search)
    }}




  #3rd Taxonomy lookup:
  DF$genus_phyletic.status<-NA
  DF$subtribe_phyletic.status<-NA
  DF$tribe_phyletic.status<-NA
  DF$subfamily_phyletic.status<-NA
  DF$family_phyletic.status<-NA
  DF$superfamily_phyletic.status<-NA
  DF$order_phyletic.status<-NA
  DF$class_phyletic.status<-NA



  ranks<-randtip::randtip_ranks()
  for (rank in ranks){
    groups<- unique(DF[,rank])
    groups<- randtip::notNA(groups)

    if (length(groups)>0){
      cat( paste0("Checking phyletic status at ", rank, " level...\n"))

      cat(paste0("0%       25%       50%       75%       100%", "\n",
                 "|---------|---------|---------|---------|", "\n"))

      for(group in groups){

        type<- randtip::MDCC.phyleticity(info, tree, MDCC.info = list("rank"= rank, "MDCC"= group))
        DF[which(DF[,rank]==group), paste0(rank,"_phyletic.status")]<-type


        v<- seq(from=0, to=40, by=40/length(groups))
        v<-ceiling(v)
        v<- diff(v)
        cat(strrep("*", times=v[which(groups==group)]))

        if(group ==groups[length(groups)]){cat("*\n")}
      }

    }

  }
  if(length(DF$Typo[DF$Typo==TRUE])>0){
    message("There may be misspelling errors in the species list or the phylogenetic tips. Please, check the outputted data frame.\n")}


  DF<-DF[,c("taxon", "PUT.status", "Typo", "Typo.names","genus", "genus_phyletic.status",
            "subtribe" , "subtribe_phyletic.status","tribe" , "tribe_phyletic.status",
            "subfamily","subfamily_phyletic.status","family","family_phyletic.status",
            "superfamily","superfamily_phyletic.status", "order","order_phyletic.status",
            "class","class_phyletic.status")]


  #4th Tree tip evaluation
  tips<- tree$tip.label

  if(length(tips[duplicated(tips)])==1){
    message("Tip ", tips[duplicated(tips)], " is duplicated in the phylogeny tips. Please remove one of them.\n")}

  if(length(tips[duplicated(tips)])>1 ){
    message("Tips ", paste0(tips[duplicated(tips)], collapse = ", "), " are duplicated in the phylogeny tips. Please remove one of them.\n")}

  subsp.tips<- tips[sapply(strsplit(tips, "_"), length)>2]

  if(length(subsp.tips)>0){
    for(ssp in subsp.tips){
      nomials<-strsplit(ssp, split="_")[[1]]
      if(paste0(nomials[1], "_", nomials[2])%in%tips &
         any(nomials[3:length(nomials)]==nomials[2])){
        message("Tips ", ssp, " and " , paste0(nomials[1], "_", nomials[2]),
                " may represent the same taxon. Please consider removing one of them.\n" )
      }
    }
  }

  #5. Tree ultrametricity evaluation
  if(!ape::is.ultrametric(tree)){
    message("The backbone tree is not ultrametric.")}


  info<- info[info$keep.tip=="1",]
  return(DF)
}



#' Convert 'info' to'input'.
#'
#' Convert an 'info' object into an 'input' one.
#'
#' @param info An 'info' data frame, including all the customized binding parameters.
#' @param tree Backbone tree.
#'
#' @return An 'input' data frame which can be fed to \code{rand.tip} alongside with a backbone tree to expand a tree.
#'
#' @author Ignacio Ramos-Gutierrez, Rafael Molina-Venegas, Herlander Lima
#'
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



