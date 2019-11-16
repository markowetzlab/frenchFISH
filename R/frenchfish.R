#get the average volume of nucleus sampled given the nucleus radius and section height
getAverageVolumeFrac<-function(r,h)
{
  Vavg=pi*h*(2/3*r^2 +r*h/3 - h^2/6)
  Vsphere=4/3*pi*r^3
  return(Vavg/Vsphere)
}

#get the maximum possible volume of nucleus sampled given the nucleus radius and section height
getMaxVolumeFrac<-function(r,h)
{
  d=0
  hh=h/2
  v=pi*hh*(r^2 -d^2 -d*hh -(1/3)*(hh^2))
  Vmax=2*v
  Vsphere=(4/3)*pi*(r^3)
  return(Vmax/Vsphere)
}

#get the minimum possible volume of nucleus sampled given the nucleus radius and section height
getMinVolumeFrac<-function(r,h)
{
  d=r-h
  Vmin=pi*h*(r^2 - d^2 - d*h - (1/3)*(h^2))
  Vsphere=(4/3)*pi*(r^3)
  return(Vmin/Vsphere)
}

#get the fraction of the nucleus sampled for a specified distance from the midpoint
getVsegFrac<-function(d,h,r)
{
  Vseg=pi*h*(r^2-d^2-h^2/12)
  Vsphere=4/3*pi*r^3
  return(Vseg/Vsphere)
}

#get volume adjusted spot count (for manual spot counting)
getCounts<-function(rawcounts,r,h)
{
  if(all(is.na(raw.counts)))
  {
    res<-NA
  }else{
    avg<-getMaxVolumeFrac(r,h)
    posterior_aqua<-MCMCpack::MCpoissongamma(rawcounts,0.01,0.01,5000)/avg
    res<-summary(posterior_aqua)$quantiles
  }
  return(res)
}

#helper function to convert spot counts and nuclear area measurements into continuous events for PP estimation
generatePPdat<-function(area,spots)
{
  indat<-0
  for(i in 1:length(area))
  {
    if(spots[i]==0)
    {
      indat[length(indat)]<-indat[length(indat)]+area[i]
    }
    else
    {
      for(k in 1:spots[i])
      {
        indat<-c(indat,area[i]/spots[i])
      }
    }
  }
  return(indat)
}

#' Main frenchFISH function for generating poisson point estimates of spot counts from
#' counts which have been automatically generated.
#' Takes a data frame where the first column is nuclear area and the remaining columns
#' are spot counts (rows are nuclear blobs).
#'
#' @param countMatrix A dataframe where the first column is nuclear area and the remaining columns are spot counts
#' @param radius The cell radius
#' @param height The section height
#' @return The \code{countMatrix} with spot counts that have been adjusted using Poisson modelling
#' @export
#' @examples
#' add(1, 1)
#' add(10, 1)
getPPcountsEstimates<-function(countMatrix,radius,height)
{
  cellarea<-pi*(radius^2)
  adjustFact<-getAverageVolumeFrac(radius,height)

  countEstimates<-c()
  for(i in 2:ncol(countMatrix))
  {
    ppindat<-generatePPdat(countMatrix$area,countMatrix[,i])
    ppcumsum<-cumsum(ppindat)/cellarea
    PPestimate<-NHPoisson::fitPP.fun(posE=ppcumsum,nobs=max(ppcumsum),start=list(b0=0),modCI=TRUE,CIty="Transf",dplot=FALSE)
    res<-data.frame(lowCI=PPestimate@LIlambda[1],median=exp(PPestimate@coef),highCI=PPestimate@UIlambda[1])
    rownames(res)<-colnames(countMatrix)[i]
    countEstimates<-rbind(countEstimates,res)
  }
  countEstimates<-countEstimates/adjustFact
  return(data.frame(countEstimates,Probe=row.names(countEstimates),row.names = NULL))
}
