#' Function to check DF0 work
#' @export
#'
#'

check.input<- function(DF0, tree, sim=0.8, verbose=F){

  DF0<- randtip::correct.DF(DF0)
  DF0[is.na(DF0$keep.tip)]<-"1"

  tree$tip.label<- gsub(" ", "_", tree$tip.label)
  DF0.taxa<-DF0$taxon
  tree.taxa<- tree$tip.label

  DF<- DF0[,c("taxon",randtip::randtip_levels())]
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

  if(length(DF$Typo[DF$Typo==TRUE])>0){
    message("There may be mistakenly written PUTs in your DF0! \nPlease check the Typo.names column in the resultant dataframe")}




  #3rd Taxonomy lookup:
  DF$genus_phyletic.status<-NA
  DF$subtribe_phyletic.status<-NA
  DF$tribe_phyletic.status<-NA
  DF$subfamily_phyletic.status<-NA
  DF$family_phyletic.status<-NA
  DF$superfamily_phyletic.status<-NA
  DF$order_phyletic.status<-NA
  DF$class_phyletic.status<-NA



levels<-randtip::randtip_levels()
  for (level in levels){
    groups<- unique(DF0[,level])
    groups<- randtip::notNA(groups)

    if(verbose){cat( paste0("Checking phyletic status at ", level, " level ... (", length(groups), " categories)"))}
    if (length(groups)>0){for(group in groups){
      type<- randtip::MDCC.phyleticity(DF0, tree, MDCC.info = list("level"= level, "MDCC"= group))
    DF[which(DF[,level]==group), paste0(level,"_phyletic.status")]<-type

      }}
   if(verbose){ cat(paste0(" Done!", "\U2713", "\n"))}
  }

DF<-DF[,c("taxon", "PUT.status", "Typo", "Typo.names","genus", "genus_phyletic.status",
          "subtribe" , "subtribe_phyletic.status","tribe" , "tribe_phyletic.status",
          "subfamily","subfamily_phyletic.status","family","family_phyletic.status",
      "superfamily","superfamily_phyletic.status", "order","order_phyletic.status",
      "class","class_phyletic.status")]


#4th Tree tip evaluation
tips<- tree$tip.label

if(length(tips[duplicated(tips)])>0){
  message("Tips ", tips[duplicated(tips)], " is duplicated in your tree!")}

subsp.tips<- tips[sapply(strsplit(tips, "_"), length)>2]

if(length(subsp.tips)>0){
  for(ssp in subsp.tips){
    nomials<-strsplit(ssp, split="_")[[1]]
    if(paste0(nomials[1], "_", nomials[2])%in%tips &
       any(nomials[3:length(nomials)]==nomials[2])){
      message("\nTips ", ssp, " and " , paste0(nomials[1], "_", nomials[2]),
              " are synonyms, and are both included in the tree" )
    }
  }
}

DF0<- DF0[DF0$keep.tip=="1",]
  return(DF)
}

#example.check<- check.input(DF0 =example , tree = tree25)

