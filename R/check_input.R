#' Function to check DF1 work
#' @export
#'
#'

check.input<- function(DF1, tree){


  DF1.taxa<-DF1$taxon
  tree.taxa<- tree$tip.label

  DF<- DF1[,c("taxon","genus", "tribe","subfamily","family", "order", "class")]
  DF$PUT.status<- NA
  DF$Name.simmilarity<- NA




  #1st we search included taxa
  for(i in 1:nrow(DF)){
    tax<- DF$taxon[i]
    if(tax %in%tree.taxa){DF$PUT.status[i]<- "Tree tip"}else{DF$PUT.status[i]<- "PUT"}
  }

  #2nd we look for name simmilarities
  for(i in 1:nrow(DF)){
    if(DF$PUT.status[i]=="PUT"){
      tax<-DF$taxon[i]
      sim.search<-tree.taxa[stringdist::stringsim(tree.taxa,tax)>0.8]
      if(length(sim.search)>0){
        sim.search<-paste0(sim.search, collapse = " / ")
        DF$Name.simmilarity[i]<- sim.search}
      rm(tax, sim.search)
    }}

  if(length(DF$Name.simmilarity[!is.na(DF$Name.simmilarity)])>0){
    warning("There may be mistakenly written PUTs in your DF1! \nPlease check the Name.simmilarity column in the resultant dataframe")}





  DF$genus_phyletic.status<-NA
  DF$tribe_phyletic.status<-NA
  DF$subfamily_phyletic.status<-NA
  DF$family_phyletic.status<-NA
  DF$order_phyletic.status<-NA
  DF$class_phyletic.status<-NA


  #3th:
levels<-c("genus", "tribe","subfamily","family", "order", "class")
  for (level in levels){
    groups<- unique(DF1[,level])
    groups<- randtip::notNA(groups)

   if (length(groups)>0){for(group in groups){
      type<- randtip::MDCC.phyleticity(DF1, tree, MDCC.info = list("level"= level, "MDCC"= group))
    DF[which(DF[,level]==group), paste0(level,"_phyletic.status")]<-type
      }}

  }

DF<-DF[,c("taxon", "PUT.status", "Name.simmilarity","genus", "genus_phyletic.status",
      "tribe" , "tribe_phyletic.status", "subfamily","subfamily_phyletic.status",
      "family","family_phyletic.status", "order","order_phyletic.status",
      "class","class_phyletic.status")]
  return(DF)
}

#example.check<- check.input(DF1 =example , tree = tree25)

