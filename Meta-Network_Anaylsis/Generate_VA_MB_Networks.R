---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Ctrl+Shift+Enter*. 

```{r}
pathS<-"P:/Surrogate/Network_Comparison/Scripts/"
source(paste0(pathS,"precision_bootstrap.R"))

pathP<-"P:/Surrogate/Network_Comparison/Data/"

freq_PSCs<- readRDS(paste0(pathP,"PSCsgllvm.RDS"))
freq_rice_gs<- readRDS(paste0(pathP,"RiceGenoSiteres_resgllvm.RDS"))
freq_DSCs<- readRDS(paste0(pathP,"DSCsgllvm.RDS"))
freq_dcsf<- readRDS(paste0(pathP,"dcsfam_.RDS"))
freq_pcsf<- readRDS(paste0(pathP,"pcsfam_gllvm.RDS"))
freq_joe<- readRDS(paste0(pathP,"rt_joe_2016gllvm.RDS"))
freq_arb2014<- readRDS(paste0(pathP,"rt_arb_2014gllvm.RDS"))
freq_arb2016<- readRDS(paste0(pathP,"rt_arb_2016gllvm.RDS"))
freq_SCGs<- readRDS(paste0(pathP,"SCGs_gllvm3.RDS"))
freq_grapes<- readRDS(paste0(pathP,"grape_newtc2.RDS"))
freq_pnas<- readRDS(paste0(pathP,"EndoRhizos_restc2.2gllvm.RDS"))

# Use for renames of PNAS data
pnas_phylos <- readRDS(file=paste0(pathP,"pnas_phylos_fam_rn.RDS"))


gllvmslist<- c(freq_SCGs,freq_pnas,freq_rice_gs, freq_grapes,
            freq_pcsf,freq_dcsf,freq_arb2014,freq_arb2016,freq_joe)
lvstouse<- lapply(gllvmslist,function(x) x$Final_gllvm$GLLVM$num.lv)


```

```{r}
freq_ricegs_pres<- lapply(freq_rice_gs,function(x)bootpres(x,bootstraps = 10000))
freq_arb_pres<- lapply(freq_arb2014,function(x)bootpres(x,bootstraps = 10000))
freq_arb2016_pres<- lapply(freq_arb2016,function(x)bootpres(x,bootstraps = 10000))
freq_joe_pres<- lapply(freq_joe,function(x)bootpres(x,bootstraps = 10000))
freq_PSCs_pres<- lapply(freq_PSCs,function(x)bootpres(x,bootstraps = 10000))
freq_DSCs_pres<- lapply(freq_DSCs,function(x)bootpres(x,bootstraps = 10000))
freq_SCGs_pres<- lapply(freq_SCGs,function(x)bootpres(x,bootstraps = 10000))

freq_grapes_pres<- lapply(freq_grapes,function(x)bootpres(x,bootstraps = 10000))
freq_pnas_pres<- lapply(freq_pnas,function(x)bootpres(x,bootstraps = 10000))

freq_dcsf_pres<- lapply(freq_pcsf,function(x)bootpres(x,bootstraps = 10000,type=2))
freq_pcsf_pres<- lapply(freq_dcsf,function(x)bootpres(x,bootstraps = 10000,type=2))

preslist<- c(freq_SCGs_pres,freq_pnas_pres,freq_ricegs_pres, freq_grapes_pres,
            freq_pcsf_pres,freq_dcsf_pres,freq_arb_pres,freq_arb2016_pres,freq_joe_pres)



saveRDS(preslist,file=paste0(pathP,"freq_pres_list.RDS"))


```


```{r}

mber<- function(themodel,type=1){ 
  tryCatch({
  if(type==1){
cov1<- themodel$cov_matrix[[1]]
  } else if(type==2){
    cov1<- themodel$cov_matrix
  }
  print(names(themodel))
lmax <- getMaxCov(cov1)
lams <- getLamPath(lmax, lmax*.001, len=100)
hugeargs <- list(lambda=lams, verbose=TRUE,lambda.min.ratio=1e-3,scr=TRUE,scr.num=(ncol(cov1)-1))
out.mb    <- pulsar(cov1, fun=huge, fargs=hugeargs, rep.num=30,
                   criterion='stars',lb.stars = TRUE,ub.stars = TRUE)
refit.mb<- refit(out.mb)
optimal<- out.mb$stars$opt.index
adj<- SpiecEasi::symBeta(refit.mb$est$beta[[optimal]])
refit.mb<- c(refit.mb,adj,optimal)
names(refit.mb)[4:5]<-c("adjacency_matrix","optimal_lambda")
return(refit.mb)},error=function(e)NA)
}

```

```{r}
freq_arb_mb<- lapply(freq_arb2014,function(x)mber(x))


freq_ricegs_mb<- lapply(freq_rice_gs,function(x)mber(x,type=2))
freq_arb2016_mb<- lapply(freq_arb2016,function(x)mber(x))
freq_joe_mb<- lapply(freq_joe2016,function(x)mber(x))
freq_PSCs_mb<- lapply(freq_PSCs,function(x)mber(x))
freq_DSCs_mb<- lapply(freq_DSCs,function(x)mber(x))
freq_SCGs_mb<- lapply(freq_SCGs,function(x)mber(x))

freq_dcsf_mb<- lapply(freq_pcsf,function(x)mber(x,type=2))
freq_pcsf_mb<- lapply(freq_dcsf,function(x)mber(x,type=2))


freq_grapes_mb<- lapply(freq_grapes,function(x)mber(x))
freq_pnas_mb<- lapply(freq_pnas,function(x)mber(x))

  #freq_DSCs_mb,

mblist<- c(freq_SCGs_mb3,freq_pnas_mb,freq_rice_gs_mb, freq_grapes_mb,
            freq_pcsf_mb,freq_dcsf_mb,freq_arb_mb,freq_arb2016_mb,freq_joe_mb)

saveRDS(mblist,file=paste0(pathP,"freq_mbs_list.RDS"))



empties<- unlist(lapply(mblist,function(x)
  if(sum(is.na(x))==0) return(sum(x$adjacency_matrix != 0))
  else return(0)))

reruns<- names(which(empties==0))


 rerunlist<- which(names(gllvmslist) %in% reruns)
rerunners<- gllvmslist[rerunlist]

#first 3 split 




```

