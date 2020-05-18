detect_negative_slope_Lowess_03 <- function(xx,yy){
require(smooth)
# N=length(yy);
#
# WindowsNumber=10;#100;
# WindowsLimits=linspace(min(xx),max(xx),WindowsNumber+1);


Smoothed_yy<- smooth::sma(yy)
Smoothed_yy<- as.vector(Smoothed_yy$fitted)
# LocalSlope=diff(yy)./diff(xx);

LocalSlope <-base::diff(Smoothed_yy)/base::diff(xx)

xc1<- which(LocalSlope > 0)
if (length(xc1)==0){ 
  xc <- xx[1] 
}else {
  xc<- xx[xc1[length(xc1)]]}
return(xc)
} 



  
