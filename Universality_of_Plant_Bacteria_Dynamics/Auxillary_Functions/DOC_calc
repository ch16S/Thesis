
func_Cal_Overlap_rJSD_from_relative_abundance <- function(physeq,NULL_COM=TRUE){
  require(phyloseq)
  # THIS FUNCTION ASSUMES A IS ALREADY NORMALIZED
if(NULL_COM==FALSE){
  if(phyloseq::taxa_are_rows(physeq)==TRUE){
    physeq<- otu_table(physeq)
    physeq<- t(physeq)
    physeq<- data.frame(physeq)
  }else{
    physeq<- otu_table(physeq)
    physeq<- data.frame(physeq)
  }
} 
  NumSamples<- nrow(physeq)
  
  ## Prepare analysis
  # preallocate memory
  
  sbj1<-sbj2<- matrix(ncol=NumSamples,nrow=NumSamples)
  for (i in 1:(NumSamples)){
    for (j in 1:NumSamples){
      sbj1[i,j]<- i
      sbj2[i,j]<- j
    }}
  
 sbj1vector<- as.vector(sbj1)
 sbj2vector<- as.vector(sbj2)

  

  ## Analize all sample pairs
  ## Make a list of sample pairs to lapply over 
  df1<- data.frame(sbj1vector,sbj2vector)
  df.list <- split(df1, seq(nrow(df1)))

thecalc<- function(i,j){ 
  # Define Kullback-Leibler Divergence
  KLD <- function(x,y) {sum(x*log(x/y))}
  # Define Jason-Shanon Divergence
  rJSD <- function(x,y){ sqrt(0.5 * KLD(x, ((x+y)/2)) + 0.5 * KLD(y, ((x+y)/2)))}  
  
SharedSpecies<-which(physeq[i,]> 0 & physeq[j,]> 0)
ifelse(length(SharedSpecies)==0,{
  Overlap <- 0
  RootJSD <- 1},{
  Overlap<- 0.5*(sum(physeq[i,SharedSpecies])+ 
                   sum(physeq[j,SharedSpecies]))
  # Renormalize the shared species
  RenormalizedA <- physeq[i,SharedSpecies]/sum(physeq[i,SharedSpecies])
  RenormalizedB <- physeq[j,SharedSpecies]/sum(physeq[j,SharedSpecies])
  # Calculate the rJSD dissimilarity
  RootJSD <- rJSD(RenormalizedA,RenormalizedB)
  })
  res<- c(Overlap,RootJSD)
  return(res)
}


results1<- lapply(df.list,function(x){
  i<- x[1,1]
  j<- x[1,2]
  thecalc(i,j)
})


Overlap<- matrix(unlist(lapply(results1,function(x) x[1])),ncol=NumSamples)
rJSD<- matrix(unlist(lapply(results1,function(x) x[2])),ncol=NumSamples)
res<- list(Overlap,rJSD)
names(res)<- c("Overlap","jsd")
return(res)}
  
 
