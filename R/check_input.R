#' Function to check info work
#' @export
#'
#'

check.info<- function(info, tree, sim=0.8, verbose=F){

  info<- randtip::correct.DF(info)
  info[is.na(info$keep.tip)]<-"1"

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

  if(length(DF$Typo[DF$Typo==TRUE])>0){
    message("There may be mistakenly written PUTs in your info! \nPlease check the Typo.names column in the resultant dataframe")}




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
      if(verbose){cat( paste0("Checking phyletic status at ", rank, " level ... (", length(groups), " categories)\n"))

      cat(paste0("0%       25%       50%       75%       100%", "\n",
                 "|---------|---------|---------|---------|", "\n"))}
      for(group in groups){
        if(group ==groups[1]&verbose){cat("*")}
      type<- randtip::MDCC.phyleticity(info, tree, MDCC.info = list("rank"= rank, "MDCC"= group))
      DF[which(DF[,rank]==group), paste0(rank,"_phyletic.status")]<-type


      if(verbose){if(length(groups)<40){
        v<- seq(from=0, to=40, by=40/length(groups))
        v<-ceiling(v)
        v<- diff(v)
        cat(strrep("*", times=v[which(groups==group)]))
      }else{
      if(which(groups==group)%in% ceiling(seq(from=0, to=length(groups),
                                              by=(length(groups)/40)))){cat("*")}}}
      }
      if(group ==groups[length(groups)]&verbose){cat("\n")}
      }

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

#5. Tree ultrametricity evaluation
if(!ape::is.ultrametric(tree)){
  message("Specified tree is not ultrametric.")}


info<- info[info$keep.tip=="1",]
  return(DF)
}

#example.check<- check.input(info =example , tree = tree25)
