
phyloseq_sep_variable <- function(physeq, variable, drop_zeroes = T){
  
  # require(phyloseq)
  # require(plyr)
  
  ## Check the input
  if(is.null(phyloseq::sample_data(physeq, errorIfNULL = F))){
    stop("Sample data is missing in the phyloseq-object.\n")
  }
  
  ## Extract samle meta-data
  mtd <- as(object = phyloseq::sample_data(physeq), Class = "data.frame")
  
  if(!variable %in% colnames(mtd)){
    stop("Grouping variable is missing from the sample data of phyloseq-object.\n")
  }
  
  if(class(mtd[, variable]) %in% c("integer", "numeric") ){
    if( length( unique(mtd[, variable]) ) > 5){
      stop("Groupping variable is numeric and it has too many levels. Consider transforming it to factor.\n")
    } else {
      warning("Groupping variable is numeric and it was coerced to factor.\n")
      mtd[, variable] <- factor(mtd[, variable])
    }
  }
  
  if(length(table(mtd[, variable])) == 1){
    cat("Warning: there is only one group of samples in the resulting list.\n")
  }
  
  ## Add sample IDs to the meta-data
  smp <- data.frame(
    SID = phyloseq::sample_names(physeq),
    mtd,
    stringsAsFactors = F)
  
  ## Exatract sample names by the specified variable
  svv <- plyr::dlply(.data = smp, .variables = variable, .fun = function(z){ z$SID })
  
  ## Extract samples by groupping variable
  res <- plyr::llply(.data = svv, .fun = function(z){ phyloseq::prune_samples(z, x = physeq) })
  
  ## Remove taxa with zero abundance
  if(drop_zeroes == TRUE){
    res <- plyr::llply(.data = res, .fun = function(x){ phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) })
  }
  
  return(res)
}


"P:\Functional_overlap2"


"H:/cbal5115/Documents/TGS/files/deML_16S/"
library(ggplot2)
library(dada2)

path <-  "P:/Functional_overlap2/"
pathF <-  "P:/Functional_overlap/filered"
pathTax <-  "P:/Functional_overlap/"

# creates a path to your files or setwd() to the location of your files (the former can create issues at times if not done correctly)
fns <- sort(list.files(path, full.names = TRUE)) # or fns <- sort(list.files( full.names = TRUE))


expectedErrors <- c(3,5)
truncationlengths <- c(270,220) 
ID<- "Bait_etal"

fnFs <- fns[grepl("_1.fastq.gz",fns)]
fnRs <- fns[grepl("_2.fastq.gz",fns)]

names<- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
namesR<- sapply(strsplit(basename(fnRs), "_"), `[`,1)
namesU <- make.unique(names)
names(fnFs) <-names
names(fnRs)<-names

namesFirst<- names

filtFs<- file.path(path,
                   "filtered", paste0(namesU, "_F_filt.fq.gz"))
filtRs<- file.path(path, 
                   "filtered", paste0(namesU, "_R_filt.fq.gz"))

out2<- filterAndTrim(fnFs,filtFs,
                    fnRs,filtRs ,
                    maxEE=c(3,5), 
                    truncLen= c(280,230),
                    compress=TRUE,verbose=TRUE,
                    rm.phix = TRUE)

saveRDS(out,paste(path,ID,"ReadsInOut.RDS",sep=""))
out3<- rbind(out,out2)

pathF<- paste0(path,"filtered/")
filts <- sort(list.files(pathF, full.names = TRUE)) # or fns <- sort(list.files( full.names = TRUE))
filtFs<- filts[grepl("_F_filt.fq.gz",filts)]
filtRs<- filts[grepl("_R_filt.fq.gz",filts)]

keep<- which(out2[,"reads.out"] > 50) # Or other cutoff

#outnames<- sapply(strsplit(basename(rownames(out3)),"_"),"[",1)
#samNames<- sapply(strsplit(basename(filtFs),"_"),"[",1)
#samNamesR<- sapply(strsplit(basename(filtRs),"_"),"[",1)
#out4<- out3[match(samNames,outnames),]
keep<- which(out2[,"reads.out"] > 50)
names<- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
namesR<- sapply(strsplit(basename(filtRs), "_"), `[`,1)

names(filtFs)<- names
names(filtRs)<- names

filtFs2<- filtFs[keep]
filtRs2<- filtRs[keep]
errF <- learnErrors(filtFs2,nbases=1e8,randomize=TRUE,multithread = TRUE)
errR <- learnErrors(filtRs2,nbases=1e8,randomize=TRUE,multithread = TRUE)
errors<- list(errF,errR)
names(errors)<-c("fwd","rvs")
saveRDS(errors,file=paste0(path,"bai_errors.RDS"))
errs<-saveRDS(file=paste0(path,"bai_errors.RDS"))
errF<- errs$fwd
errR<- errs$rvs
mergers<- vector("list", length(filtFs2)) 
names(mergers) <- names(filtFs2)
for(sam in 1:length(filtFs2)) {
cat("Processing:", sam, "\n")
  tryCatch({
  derepF <- derepFastq(filtFs2[[sam]])
  ddFs <- dada(derepF, err=errF, multithread=FALSE)
  derepR <- derepFastq(filtRs2[[sam]])
  ddRs <- dada(derepR, err=errR, multithread=FALSE)
  merger<- mergePairs(ddFs,derepF,ddRs,derepR )},error=function(e)NA)
  mergers[[sam]] <- merger
}

remove<- which(is.na(names(mergers)))

seqtab<-makeSequenceTable(mergers[-remove])

saveRDS(seqtab,paste(path,ID,"seqtab.RDS",sep=""))
no.chim<- removeBimeraDenovo(seqtab)

taxgg<- assignTaxonomy(colnames(no.chim),paste0(pathTax,"gg_13_8_train_set_97.fa.gz"))

meta<- read.csv(paste0(path,"bai_etal_meta.csv"),header=TRUE,row.names = 1)
meta2<- meta[,1:5]

startingtax<- read.csv("P:/Surrogate/bai2015_taxa_used.csv",header=TRUE,row.names=NULL)

saveRDS(no.chim,paste(path,ID,"noChim.RDS",sep=""))
#

ps<- phyloseq(otu_table(no.chim,taxa_are_rows=FALSE),tax_table(taxgg))

saveRDS(ps,paste(path,ID,"ps.RDS",sep=""))

filterPhyla = c("f__mitochondria")
surr1mito<- subset_taxa(ps,Family %in% filterPhyla)
filterPhyla = c("c__Chloroplast")
surr1chloro<-  subset_taxa(ps,Class %in% filterPhyla)
surr1organeller<- rowSums(otu_table(surr1mito))+rowSums(otu_table(surr1chloro))

ps2<- subset_taxa(ps, ! Class %in% c("c__Chloroplast"))
ps3<- subset_taxa(ps2, ! Family %in% c("f__mitochondria"))
ps4<- subset_taxa(ps3,!is.na(Phylum))

plastidRatio<- log((sample_sums(ps4)+1)/(surr1organeller+1))
ps5<- filter_taxa(ps4,function(x)sum(x>2)>1,TRUE)
ps6<- subset_samples(ps5_2,!Compartment %in% c("Leaf","Control"))
psroot<- subset_samples(ps5_2,Compartment %in% c("Root"))

bytreat<- phyloseq_sep_variable(ps6,"Treatment.1")

keep<- which(unlist(lapply(bytreat,function(x)length(unique(sample_data(x)$Compartment))))>1)

theRoots<- bytreat[keep]
theRoots2<- lapply(theRoots,function(x)prune_samples(sample_sums(x)>500,x))


theRoots2<- lapply(theRoots2,function(x)tax_glom(x,"Genus"))
theDonors<- lapply(theRoots2,function(x)prune_samples(sample_data(x)$Compartment =="Ino",x))
theDonorstest<- lapply(theDonors,function(x)prune_taxa(Family %in% startingtax$gg_Genus,x))

results<- vector("list",length=length(theRoots))
for(i in 1:length(theRoots)){
  donor<- subset_samples(theRoots[[i]],Compartment %in% "Ino" )
  hosts<- subset_samples(theRoots[[i]],!Compartment %in% "Ino" )
  donor2<- filter_taxa(donor,function(x)sum(x>0)>0,TRUE)
  donor3<- subset_taxa(donor2,Family %in% startingtax$gg_family)
  target<- subset_taxa(hosts, taxa_names(hosts) %in% taxa_names(donor3))
  nontarget<- subset_taxa(hosts,! taxa_names(hosts) %in% taxa_names(donor3))
  targetsum<- sample_sums(target)
  nontargetsum<- sample_sums(nontarget)
  logratio<- log((targetsum+1)/(nontargetsum+1))
  targetr<- estimate_richness(target,measures="Observed")
  nontargetr<- estimate_richness(nontarget,measures="Observed")
  logRratio<- log((targetr$Observed+1)/(nontargetr$Observed+1))
  donormerge<- merge_samples(donor3,"Compartment" , function(x)sum(x))
  donorr<- estimate_richness(donormerge,measures="Observed")
  recovered<- targetr$Observed/donorr$Observed
  df1<- data.frame(target=targetsum,nontarget=nontargetsum,
                   DR_ratio=logratio,target_r=targetr$Observed,
                   nontarget_r=nontargetr$Observed,RR_ratio=  logRratio,
                   spp_recoved=recovered,donor_richness=donorr$Observed)
  results[[i]]<- df1
}
names(results)<- names(theRoots2)
final_df<- plyr::ldply(results,as.data.frame,by=.id)
final_df

saveRDS(final_df,paste0(pathContam,"Bai_etal_transplantres.RDS"))





