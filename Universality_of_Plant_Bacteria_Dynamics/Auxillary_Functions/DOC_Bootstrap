

func_cal_rlowess_bootstrap <- function(x ,NumBootstraps= 200,span=0.2){
  require(phyloseq)
  require(pracma)
  Overlap<- x$Overlap
  RootJSD<- x$jsd 
  NumSamples <- ncol(Overlap)
  diag(Overlap)<- NA
  diag(RootJSD)<- NA
  

  sbj1<-sbj2<- matrix(ncol=NumSamples,nrow=NumSamples)
  for (i in 1:(NumSamples)){
    for (j in 1:NumSamples){
      sbj1[i,j]<- i
      sbj2[i,j]<- j
    }}
  
  xs <- as.vector(na.omit(as.dist(Overlap)))
 # xs<- na.omit(xs)
  # [xs,ind]=sort(xs);
  if(NumSamples < 50){
    xs <- seq(from=min(xs),to=max(xs),length.out=20)
    
  }else{
  xs <- seq(from=min(xs),to=max(xs),length.out=100)}
  
  
  p_lmes <- vector(length=NumBootstraps)
  y_Bootstrap <- matrix(ncol=length(xs),nrow=NumBootstraps)
  for(i in 1:NumBootstraps){
    # resample NumSamples samples
    k <- sample.int(NumSamples,NumSamples,replace=TRUE)
    ChosenOverlap <- Overlap[k,k]##ok<*PFBNS>
    ChosenRJSD <- RootJSD[k,k]
    chosn_sbj1<- sbj1[k,k]
    chosn_sbj2<- sbj2[k,k]
    OverlapVector <-   as.vector(na.omit(as.dist(ChosenOverlap)))
    RJSDVector <-  as.vector(na.omit(as.dist(ChosenRJSD)))
    chosn_sbj1v <- as.vector(na.omit(as.dist( chosn_sbj1)))
    chosn_sbj2v <-   as.vector(na.omit(as.dist(chosn_sbj2)))
    
    
    #     ys=mylowess([OverlapVector  RJSDVector],xs,span);
    xy1<- matrix(c(OverlapVector,  RJSDVector),ncol=2)
    if(sd(xy1[,1])==0){
      ys<- rep(1,length(xs))
    }else{
      ys <- myrlowess(xy1,xs,span)
    }
    
    y_Bootstrap[i,] <- ys
    
    
    tokeep<- which( OverlapVector > median(as.vector(na.omit(as.dist(Overlap)))))
    chosn_sbj1v <- as.factor(chosn_sbj1v[tokeep])
    chosn_sbj2v <-  chosn_sbj2v[tokeep]
    RJSDVector<- RJSDVector[tokeep]
    OverlapVector<-  OverlapVector[tokeep]
    
    
    
    lme<- lme4::lmer( RJSDVector~ 1 + OverlapVector +(1|chosn_sbj1v)+(1|chosn_sbj2v))
    
    p_lmes[i]<-summary(lme)$coefficients[2]
  }# for i
  
  bootstrap_res<- list(y_Bootstrap,xs,p_lmes)
  names(bootstrap_res)<- c("Y-Bootstrap","xs","neg_slopes")
  return(bootstrap_res)
  
}# func_cal_rlowess_bootstrap


######################
myrlowess <- function(xy,xs,span=0.3){
  require(pracma)
  #MYLOWESS Lowess smoothing, preserving x values
  #   YS=MYLOWESS(XY,XS) returns the smoothed version of the x/y data in the
  #   two-column matrix XY, but evaluates the smooth at XS and returns the
  #   smoothed values in YS.  Any values outside the range of XY are taken to
  #   be equal to the closest values.
  
  
  # Sort and get smoothed version of xy data
  xy<- xy[order(xy[,1],decreasing = FALSE ),]
  
  x1 <- xy[,1]
  y1 <- xy[,2]
  ys1 <- lowess(x1,y1,f=span)
  
  # Remove repeats so we can interpolate
  removal<- NULL
  x_diff<- diff(x1)
  for(i in 1:length(x_diff)){
    if(x_diff[i]==0){
      t<- i
      removal<- c(removal,t)
    }
  }
  
  if(length(removal)> 0){
  removal<- removal+1 
  x1<- x1[-removal]
  ys2<- ys1$y[-removal]
  }else{
    ys2<- ys1$y
  }
  
  
  # Interpolate to evaluate this at the xs values
  
  ys <-  approx(x1,ys2,xout=xs,method="linear")
  #ys<- pracma::interp1(x1,ys2,xs,method="linear")
  ys3<- ys$y
  #ys3<- ys
  # Some of the original points may have x values outside the range of the
  # resampled data.  Those are now NaN because we could not interpolate them.
  # Replace NaN by the closest smoothed value.  This amounts to extending the
  # smooth curve using a horizontal line.
  
  undershoot<- which(xs < x1[1])
  firstval<- undershoot[length(undershoot)]+1
  ys3[undershoot]<- ys3[firstval]


overshoot<- which(xs > (x1[length(ys2)]))
lastval<- overshoot[1]-1
ys3[overshoot]<- ys3[lastval]

return(ys3)

}# function myrlowess


