
# Read in the data, remove plant organellar DNA 

library(phyloseq);library(gridExtra);library(ggplot2);library(plyr);library(dplyr)
pathPlots<- "P:/Surrogate/Transplant_Anaylsis/Plots/"

path<- "P:/USYD_Stuff/Surrogate/files/"
asv<- readRDS( paste0(path,"surr_ASV_V3.RDS"))
tax<- readRDS("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/surr_TAX_V3.RDS")
meta<-read.csv("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/surr_metav2.csv",header = TRUE,row.names=1)
tree<- readRDS("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/fitGTR.RDS")

path<- "P:/USYD_Stuff/Surrogate/files/"

ps<-phyloseq(otu_table(asv,taxa_are_rows=FALSE),
             tax_table(tax),sample_data(meta),phy_tree(tree$tree))

filterPhyla = c("f__mitochondria")
ps0<- subset_taxa(ps,!Family %in% filterPhyla)
filterPhyla = c("c__Chloroplast")
ps0<- subset_taxa(ps0,!Class %in% filterPhyla)
ps0 <- subset_taxa(ps0, !is.na(Phylum) & !Phylum %in% c("", " uncharacterized"))

table(tax_table(ps0)[, "Family"], exclude = NULL) 
table(tax_table(ps0)[, "Class"], exclude = NULL) 

remove<-which(rowSums(otu_table(ps0))<500)
toremove<-c(names(remove))
ps0<-subset_samples(ps0,!rownames(otu_table(ps0)) %in% toremove )

surr1<-subset_samples(ps0,sample_data(ps0)$Experiment=="Surrogate")

asvs<-readRDS("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/surr2v3_asv.RDS") ## 
tax<-readRDS("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/SurrRun2_tax.RDS") ## 
samp<- read.csv("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/surr2_sData.csv",header=TRUE,row.names=1) 

ps2 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", " uncharacterized"))
mito<-c("f__mitochondria")
ps2<- subset_taxa(ps2,!Family %in% mito)
chloro = c("c__Chloroplast")
surr2<- subset_taxa(ps2,!Class %in% chloro) ## note the change in taxonomic level to remove an order opposed to family (applies to any level)

surr_merge<- merge_phyloseq(surr1,surr2)



#### filter function 

boral_filter<- function(physeq,howprev=.1,howabund=0.01){
  prev<- nsamples(physeq)*howprev
  abund<- (summary(rowSums(otu_table(physeq)))[[1]]*howabund)
  physeq1<- filter_taxa(physeq,function(y)sum(y>0)>prev,TRUE)
  physeq2<- filter_taxa(physeq1,function(z)sum(z>abund)>0,TRUE)
  return(physeq2)
}


#### merge ASVs from different runs, subset data

surr_merge<- subset_samples(surr_merge,sample_data(surr_merge)$Donor.Surrogate=="Surrogate")



rootsasv<- subset_samples(surr_merge,Compartment %in% "Root")
leavesasv<- subset_samples(surr_merge,Compartment %in% "Leaf")
rootasvf <-boral_filter(rootsasv)
leafasvf <-boral_filter(leavesasv)


library(boral)
library(phyloseq)
library(R2jags)
library(plyr)

## Function to make model matrix 

deviate<- function(x,des=TRUE){
  sdata<- sample_data(x)
  sdata<- data.frame(sdata)
  sdata$Time<- as.factor(sdata$Time)
  contrasts(sdata$Time)= contr.sum
  contrasts(sdata$Plant)=contr.sum
  contrasts(sdata$Donor_Species)= contr.sum
  contrasts(sdata$Cultured.Total)= contr.sum
  deviate<- model.matrix(~Plant*Donor_Species*Cultured.Total*Time,data=sdata)
  if(des==TRUE){ # no intercept
    covarsdev<- deviate[,c(-1)]
  }else{ # intercept 
    covarsdev<- deviate
    
  }
  return(covarsdev)}



rdevcovars<-deviate(rootsasvf)
ldevcovars<-deviate(leafasvf)

mcmcs= list(n.burnin=50000, n.iteration=210000,n.thin=10,seed=123)
priors = list(type = c("normal","normal","normal","uniform"),hypparams = c(10, 10,6.25, 30),ssvs.index = -1)


rtcr<- boral(otu_table(rootasvf),X=rdevcovars,family="negative.binomial",num.lv=3 ,row.eff="fixed",
             row.ids=matrix(1:nsamples(roots),ncol=1),mcmc.control = mcmcs,
             prior.control = priors, save.model=TRUE)

saveRDS(rtplv ,file="boral_root_results/rt_asv_cr_surr.RDS")


rtplv<- boral(otu_table(rootasvf),family="negative.binomial",num.lv=3 ,row.eff="fixed",
             row.ids=matrix(1:nsamples(roots),ncol=1),mcmc.control = mcmcs,
             prior.control = priors, save.model=TRUE)

saveRDS(rtplv ,file="boral_root_results/rt_asv_cr_plv.RDS")

lfcr<- boral(otu_table(leafasvf ),X=ldevcovars,family="negative.binomial",num.lv=3 ,row.eff="fixed",
             row.ids=matrix(1:nsamples(leafasvf ),ncol=1),mcmc.control = mcmcs,
             prior.control = priors, save.model=TRUE)

saveRDS(lfcr,file="boral_root_results/lf_asv_cr_surr.RDS")


lfplv<- boral(otu_table(leafasvf ),family="negative.binomial",num.lv=3 ,row.eff="fixed",
             row.ids=matrix(1:nsamples(leafasvf ),ncol=1),mcmc.control = mcmcs,
             prior.control = priors, save.model=TRUE)

saveRDS(lfplv ,file="boral_root_results/lf_asv_plv_surr.RDS")



par(mfrow=c(4,4))
plot("DS_Residuals_gllvms.png")
plot(rtcr)
plot(rtplv)
plot(lfcr)
plot(lfplv)

dev.off()

plot("Convergence_gllvms.png")

mcmcs1<- get.mcmcsamples(rtcr)
toplot<- sample(1:ncol(mcmcs1), ncol(mcmc)*.1)
plot(mcmcs1[,toplot])

mcmcs1<- get.mcmcsamples(rtplv)
toplot<- sample(1:ncol(mcmcs1), ncol(mcmc)*.1)
plot(mcmcs1[,toplot])

mcmcs1<- get.mcmcsamples(lfcr)
toplot<- sample(1:ncol(mcmcs1), ncol(mcmc)*.1)
plot(mcmcs1[,toplot])

mcmcs1<- get.mcmcsamples(lfplv)
toplot<- sample(1:ncol(mcmcs1), ncol(mcmc)*.1)
plot(mcmcs1[,toplot])

dev.off()
