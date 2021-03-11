#' Function to create DF1 given a species vector(column)
#' @export

build.input<- function(species,find.MDCC=TRUE, MDCC.rank="family",db="ncbi", mode="list",tree){


  if(is.data.frame(species)){
  if(ncol(species)!=1){stop("Species list must be a vector or a single-column dataframe!")}else{species<- species[,1]}}
  if(!(is.vector(species) )){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.list(species) ){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.null(tree)){stop("A tree must be provided.")}
  if (!inherits(tree, "phylo")) { stop("tree should be an object of class \"phylo\".")}


  names_df <- c("taxon", "genus", "phyleticity",  "MDCC", "other.MDCC", "aggregate.subspecies",
                "relative.species","synonim.genus","sibling.genus", "forbidden.branches")

  DF1<- as.data.frame(matrix(nrow = length(species), ncol = length(names_df)))
  names(DF1)<- names_df




  species<- gsub(" ", "_", species)

  DF1$taxon<- species

  DF1$genus<- word(species, 1, sep="_")



  if(mode=="phylomatic"){

    DF2<- as.data.frame(matrix(nrow = length(tree$tip.label), ncol = 3))
    names(DF2)<- c("taxon","genus", "MDCC")

    DF2$taxon<- tree$tip.label

    DF2$genus<- word(tree$tip.label, 1, sep="_")
  }


  if(isTRUE(find.MDCC)){
    if(mode=="phylomatic"){genera<- unique(c(word(species, 1, sep="_"), word(tree$tip.label, 1, sep="_")))}else{
      genera<- unique(word(species, 1, sep="_"))}

    for(i in 1:length(genera)){
      tryCatch({
        search <- taxize::classification(as.character(genera[i]), db = db)[[1]]
        MDCC<-search[which(search$rank==MDCC.rank), "name"]
        DF1$MDCC[DF1$genus==genera[i]]<- MDCC
        if(mode=="phylomatic"){DF2$MDCC[DF2$genus==genera[i]]<- MDCC}
      }, error=function(e){
        # Assign NA to fetching errors
        DF1[DF1$genus==genera[i], "MDCC"] <- NA
        if(mode=="phylomatic"){DF2[DF2$genus==genera[i], "MDCC"] <- NA}
      })

      # Avoid ip blocks. Taxize allows only 3 searches per second.
      Sys.sleep(0.33)
    }}


  #search for genera phyletic statuses

  genera<- unique(DF1$genus)
  for(i in 1:length(genera)){
    genus<- genera[i]
    phyleticity<-randtip::phyleticity(tree = tree, genus=genus)
    DF1[DF1$genus==genus,]$phyleticity<- phyleticity
    rm(genus, phyleticity)
  }


  if(mode=="phylomatic"){return(list("DF1"=DF1, "DF2"=DF2))}else{return(DF1)}
}

#example<-build.input(species = phylo25.table$taxon, find.MDCC = T , mode = "list", tree=tree25)
