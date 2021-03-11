#' Function to check DF1 work
#' @export
#'
#'

check.input<- function(DF1, tree){


  DF1.taxa<-DF1$taxon
  tree.taxa<- tree$tip.label

  DF<- DF1[,c("taxon","genus", "phyleticity","MDCC","other.MDCC")]
  DF$PUT.status<- NA
  DF$Name.simmilarity<- NA
  DF$Using.MDCC<- NA
  DF$Level.MDCC<- NA



  #1st we search included taxa
  for(i in 1:nrow(DF)){
    tax<- DF$taxon[i]
    if(tax %in%tree.taxa){DF$PUT.status[i]<- "Tree tip"}else{DF$PUT.status[i]<- "PUT"}
  }

  #2nd we look for name simmilarities
  for(i in 1:nrow(DF)){
    if(DF$PUT.status[i]=="PUT"){
      tax<-DF$taxon[i]
      sim.search<-tree.taxa[stringsim(tree.taxa,tax)>0.8]
      if(length(sim.search)>0){
        sim.search<-paste0(sim.search, collapse = " / ")
        DF$Name.simmilarity[i]<- sim.search}
      rm(tax, sim.search)
    }}

  if(length(DF$Name.simmilarity[!is.na(DF$Name.simmilarity)])>0){
    warning("There may be mistakenly written PUTs in your DF1! \nPlease check the Name.simmilarity column in the resultant dataframe")}






  #3th:
  genera<- unique(DF1$genus)
  for(i in 1:length(genera)){
    genus<- genera[i]
    if(unique(DF$phyleticity[DF$genus==genus])!="Not included"){
      DF$Level.MDCC[DF$genus==genus]<-"genus"
      DF$Using.MDCC[DF$genus==genus]<- genus
    }}


  for(i in which(is.na(DF$Using.MDCC)) ){

    other.MDCC<- DF$other.MDCC[i]
    commonMDCC<- DF$MDCC[i]

    if(!is.na(other.MDCC) & length(DF[DF$other.MDCC==other.MDCC,])>1){
      DF$Level.MDCC[i]<-"other.MDCC"
      DF$Using.MDCC[i]<- other.MDCC
    }else{
      if(is.na(commonMDCC)){stop("At least one MDCC must be defined for genus ", DF$genus[i])}
      DF$Level.MDCC[i]<-"MDCC"
      DF$Using.MDCC[i]<- commonMDCC
    }
  }


  #DF[which(DF$PUT.status!="PUT"),"Level.MDCC"]<- "-"
  #DF[which(DF$PUT.status!="PUT"),"Using.MDCC"]<-"-"




  return(DF)
}

#example.check<- check.input(DF1 =example , tree = tree25)

