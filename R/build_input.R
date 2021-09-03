#' Function to create DF1 given a species vector(column)
#' @export

build.input<- function(species, tree, find.MDCC=FALSE, db="ncbi", mode="list",  genus=F){


  if(is.data.frame(species)){
  if(ncol(species)!=1){stop("Species list must be a vector or a single-column dataframe!")}else{species<- species[,1]}}
  species<- as.vector(species)
  if(!(is.vector(species) )){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.list(species) ){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.null(tree)){stop("A tree must be provided.")}
  if (!inherits(tree, "phylo")) { stop("tree should be an object of class \"phylo\".")}
  if(!(mode %in% c("list", "phylomatic"))) {stop("type must be \"list\" or \"phylomatic\" ")}



  names_df <- c("taxon",  "taxon1", "taxon2",
                randtip::randtip_levels(),
                "agg.ssp","rand.type", "poly.ins",
                "resp.mono", "resp.para", "resp.sing", "keep.tip")

  species<- gsub(" ", "_", species)

  spp.in.tree<- tree$tip.label
  spp.original<- species

  if(mode=="phylomatic"){
        species<- c(species, spp.in.tree[!(spp.in.tree%in%species)])
  }

  DF0<- as.data.frame(matrix(nrow = length(species), ncol = length(names_df)))
  names(DF0)<- names_df



  DF0$taxon<- species
  onlygenus<-which(sapply(strsplit(species, "_"), length)==1)
  if(length(onlygenus)>0){
    for(t in onlygenus){
      if(paste0(DF0$taxon[t], "_sp.")%in% tree$tip.label){
        DF0$taxon[t]<-paste0(DF0$taxon[t], "_sp2.")}else{DF0$taxon[t]<-paste0(DF0$taxon[t], "_sp.")}

    }
  }


  if(isTRUE(genus)){
    if(any(sapply(strsplit(species, "_"), length)>1)){stop("Taxa must specify only genera for \"genus\" mode")}
    if(any(sapply(strsplit(tree$tip.label, "_"), length)>1)){stop("Tree tips must represent only genera for \"genus\" mode")}
    DF0$taxon<- randtip::firstword(DF0$taxon)
  }

  DF0$genus<- randtip::firstword(species)


  genera<- unique(randtip::firstword(species))

  searching.categories<- randtip::randtip_levels()[-1]

    if(find.MDCC){for(i in 1:length(genera)){
      tryCatch({
        search <- suppressMessages(taxize::classification(as.character(genera[i]), db = db))[[1]]

        for(cat in searching.categories){
          if(length(search[which(search$rank==cat), "name"])==0){DF0[DF0$genus==genera[i], cat]<-NA}else{
            cats<-search[which(search$rank==cat), "name"]
            if(length(cats)>1){cats<-cats[1]}
            DF0[DF0$genus==genera[i], cat]<- cats}
        }


      }, error=function(e){
        # Assign NA to fetching errors
        DF0[DF0$genus==genera[i], searching.categories] <- NA
      })

      # Avoid ip blocks. Taxize allows only 3 searches per second.
      Sys.sleep(0.33)
    }}

  DF0[!(species%in%spp.original),
      c("agg.ssp","rand.type", "poly.ins", "resp.mono", "resp.para", "resp.sing" )]<-"-"
  DF0$keep.tip[!(species%in%spp.original)]<- 0
  DF0$keep.tip[  species%in%spp.original ]<- 1
  return(DF0)
}

#example<-build.input(species = phylo25.table$taxon, find.MDCC = T , mode = "list", tree=tree25)

complete.input<- function(DF0, tree, verbose=F){

  DF1<-DF0
  DF1[is.na(DF1$keep.tip)]<-"1"
  tree$tip.label <- gsub(" ", "_", tree$tip.label)
  DF1<- randtip::correct.DF(DF1)
  DF1$taxon <- gsub(" ", "_", DF1$taxon)

  tree$tip.label <- gsub("_x_", "_x-", tree$tip.label)
  tree$tip.label <- gsub("_X_", "_x-", tree$tip.label)


  DF1_search<- randtip::usingMDCCfinder(DF1 = DF1, taxon = DF1$taxon, tree = tree, verbose)

  DF1$using.MDCC     <- DF1_search[[1]]
  DF1$using.MDCC.lev <- DF1_search[[2]]
  DF1$using.MDCC.phylstat <- DF1_search[[3]]
  #3.2 Taxa with no MDCC

  not.included<- DF1[is.na(DF1$using.MDCC),]
  if(length(not.included$taxon) > 0){
    message("The following taxa do not have a defined MDCC and cannot be randomized:\n",
            paste0(not.included$taxon, "\n"))}

  return(DF1)
}
