#' Function to create info given a species vector(column)
#' @export

build.info<- function(species, tree=NULL, find.ranks=TRUE, db="ncbi", mode="backbone",
                      interactive=F,genus=F){


  if(is.data.frame(species)){
  if(ncol(species)!=1){stop("Species list must be a vector or a single-column dataframe!")}else{species<- species[,1]}}
  species<- as.vector(species)
  if(!(is.vector(species) )){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.list(species) ){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.null(tree) & mode=="backbone"){stop("A tree must be provided for \"backbone\" mode.")}
  if (!is.null(tree) & !inherits(tree, "phylo")) { stop("tree should be an object of class \"phylo\".")}
  if(!(mode %in% c("list", "backbone"))) {stop("type must be \"list\" or \"backbone\" ")}



  names_df <- c("taxon", randtip::randtip_ranks(),
               "taxon1", "taxon2","rand.type","polyphyly.scheme", "use.paraphyletic",
               "use.singleton","respect.mono", "respect.para","clump.PUTs",
               "keep.tip")

  species<- gsub(" ", "_", species)

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
        search <- suppressMessages(taxize::classification(as.character(genera[i]),
                                                          db = db, rows=Inf))[[1]]}

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
    }}

  info[!(species%in%spp.original),
      c("rand.type", "polyphyly.scheme","use.paraphyletic", "use.singleton",
        "respect.mono","respect.para","clump.PUTs" )]<-"-"
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


  input_search<- randtip::usingMDCCfinder(input = input, taxon = input$taxon, tree = tree)

  input$using.MDCC     <- input_search[[1]]
  input$using.MDCC.lev <- input_search[[2]]
  input$using.MDCC.phylstat <- input_search[[3]]
  #3.2 Taxa with no MDCC

  not.included<- input[is.na(input$using.MDCC),]
  if(length(not.included$taxon) > 0){
    message("The following taxa do not have a defined MDCC and cannot be randomized:\n",
            paste0(not.included$taxon, "\n"))}

  return(input)
}


