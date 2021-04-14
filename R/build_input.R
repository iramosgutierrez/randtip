#' Function to create DF1 given a species vector(column)
#' @export

build.input<- function(species, tree, find.MDCC=TRUE, db="ncbi", mode="list"){


  if(is.data.frame(species)){
  if(ncol(species)!=1){stop("Species list must be a vector or a single-column dataframe!")}else{species<- species[,1]}}
  if(!(is.vector(species) )){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.list(species) ){stop("Species list must be a vector or a single-column dataframe!")}
  if(is.null(tree)){stop("A tree must be provided.")}
  if (!inherits(tree, "phylo")) { stop("tree should be an object of class \"phylo\".")}


  names_df <- c("taxon", "genus", "tribe", "subfamily", "family", "order", "class",
                "agg.ssp","rand.type", "poly.ins")

  DF1<- as.data.frame(matrix(nrow = length(species), ncol = length(names_df)))
  names(DF1)<- names_df




  species<- gsub(" ", "_", species)

  DF1$taxon<- species

  DF1$genus<- randtip::firstword(species)



  if(mode=="phylomatic"){

    DF2<- as.data.frame(matrix(nrow = length(tree$tip.label), ncol = 7))
    names(DF2)<- c("taxon", "genus", "tribe", "subfamily", "family", "order", "class")

    DF2$taxon<- tree$tip.label

    DF2$genus<- stringr::word(tree$tip.label, 1, sep="_")
  }


  if(mode=="phylomatic"){genera<- unique(c(randtip::firstword(species), randtip::firstword(tree$tip.label)))}else{
      genera<- unique(randtip::firstword(species))}

  searching.categories<- c("tribe", "subfamily", "family", "order", "class")

    for(i in 1:length(genera)){
      tryCatch({
        search <- suppressMessages(taxize::classification(as.character(genera[i]), db = db))[[1]]

        for(cat in searching.categories){
          if(length(search[which(search$rank==cat), "name"])==0){DF1[DF1$genus==genera[i], cat]<-NA}else{
          DF1[DF1$genus==genera[i], cat]<- search[which(search$rank==cat), "name"]}
        }

        if(mode=="phylomatic"){
          for(cat in searching.categories){
            if(length(search[which(search$rank==cat), "name"])==0){DF2[DF2$genus==genera[i], cat]<-NA}else{
              DF2[DF2$genus==genera[i], cat]<- search[which(search$rank==cat), "name"]}
          }}


      }, error=function(e){
        # Assign NA to fetching errors
        DF1[DF1$genus==genera[i], searching.categories] <- NA
        if(mode=="phylomatic"){DF2[DF2$genus==genera[i], searching.categories] <- NA}
      })

      # Avoid ip blocks. Taxize allows only 3 searches per second.
      Sys.sleep(0.33)
    }


  if(mode=="phylomatic"){return(list("DF1"=DF1, "DF2"=DF2))}else{return(DF1)}
}

#example<-build.input(species = phylo25.table$taxon, find.MDCC = T , mode = "list", tree=tree25)
