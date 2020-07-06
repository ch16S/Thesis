
setwd("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/deML_5/")

set.seed(123)
path <- "H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/deML_5/"  # creates a path to your files or setwd() to the location of your files (the former can create issues at times if not done correctly)
fns <- sort(list.files(path, full.names = TRUE)) # or fns <- sort(list.files( full.names = TRUE))
fnFs <- fns[grepl("r1.fq.gz", fns)] ## If single end reads you do not need both fnFs and fnRs 
fnRs <- fns[grepl("r2.fq.gz", fns)]## please note r1 and r2 maybe denoted as R1 or R2 on your files and/or fastq/fq- check this

sam.names <- sapply(strsplit(fnFs, "_"), `[`, 5)
sam.names2 <- sapply(strsplit(fnFs, "_"), `[`, 6)
sampleNames<- paste(sam.names,"_",sam.names2,sep="")

dir.create("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/Illumina_2_filtered" )

pathfiltered<- "H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/Illumina_2_filtered"


filtFs <- file.path(pathfiltered, "filtered", paste0(sampleNames, "_F_filt.fq.gz"))
filtRs <- file.path(pathfiltered, "filtered", paste0(sampleNames, "_R_filt.fq.gz"))



if(length(fnFs)!=length(fnRs))stop("FWD and RVS don't match go back and work out why")                                           

## Then filter -> change parameters to fit your data. namely the length to truncate the sequences 

out<- filterAndTrim(fnFs,filtFs,
              fnRs,filtRs,maxEE=c(2,4), 
              trimLeft=c(21,17),
              truncLen = c(280,260),
              compress=TRUE,verbose=TRUE,
              rm.phix = TRUE)


sample.names<- sapply(strsplit(basename(filtFs),"_"),"[",1)
sample.namesR<- sapply(strsplit(basename(filtRs),"_"),"[",1)

if(!identical(sample.names,sample.namesR)) stop("FWD and RVS don't match")
if(!identical(sample.names,sampleNames)) stop("FWD and RVS don't match")

names(filtFs )<-sampleNames
names(filtRs )<- sampleNames



# nreads is the number of reads you use to learn errors, ~1000000 is recommended for a hiseq run, so you can halve it if you like, i.e 5e5 . 

errF <- learnErrors(filtFs , nbases=1e8, multithread=FALSE,randomize=TRUE)
errR <- learnErrors(filtRs , nbases=1e8, multithread=FALSE,randomize=TRUE)

plotErrors(errF,nominalQ = TRUE) ## Check error,, change to errR for reverse read error rates


# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sampleNames))
names(mergers) <- sampleNames
for(sam in sampleNames) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=FALSE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=FALSE)
  merger <- mergePairs(ddF, derepF, ddR, derepR) ####Add this for reads that dont overlap: justConcatenate=TRUE
  mergers[[sam]] <- merger
}

seqtab <- makeSequenceTable(mergers)

hist(nchar(getSequences(seqtab)), main="Distribution of sequence lengths")
seqtab2<- seqtab[,nchar(colnames(seqtab)) %in% seq(393,431)] ## Excise amplicons of specific size. 
saveRDS(seqtab, "H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/SurrRun2_seqtabV1.RDS") ## change to directory you wish to save 


nochim<- removeBimeraDenovo(seqtab2)

table(nchar(getSequences(nochim)))

100-((sum(nochim)/sum(seqtab2))*100) # Number of reads that are chimeric

100-(length(nochim[1,])/length(seqtab2[1,])*100)
## Should be a large number of SVs are chimeras, but a small number of reads

saveRDS(nochim, "H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/SurrRun2_noChimV1.RDS") ## change to directory you wish to save 


ref_fasta <- "H:/cbal5115/Documents/PhD/Protocols/gg_13_8_train_set_97.fa.gz" ## Use appropriate taxonomy here 
taxtab <- assignTaxonomy(nochim, refFasta = ref_fasta,tryRC=TRUE)
```
sData<- read.csv("H:/cbal5115/Documents/PhD/Probiotics/Surrogate/files/surr2sData.csv",header=TRUE,row.names=1)

###

ps<- phyloseq(otu_table(nochim2,taxa_are_rows = FALSE),tax_table(taxtab),sample_data(sData),phy_tree(fitGTR$tree))

Relabelled the samples - file nochim.csv holds the key 1
nochim2<- nochim 
Old:  
orig_names<- c("BNLFCU1" ,"QALFCU1","QALF2.3", "GABLFP4")
new_names<- c("QALF2.3", "GABLFP4","BNLFCU1" ,"QALFCU1")
nochim2<- read.csv("nochim.csv",header = TRUE,row.names = 1)

Old
BNLFCU1 
QALFCU1 
QALF2.3 
GABLFP4

New: 
QALF2.3
GABLFP4
BNLFCU1
QALFCU1

for( i in seq_along(orig_names)){
  name<- orig_names[[i]]
  newname<- new_names[[i]]
  namesMaster <-gsub(name,newname,namesMaster)
  
  
}

o__Chlorophyta


