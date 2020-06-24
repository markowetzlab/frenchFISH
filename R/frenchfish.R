#' Helper function to check if all values in the input count
#' matrix are either NA, NaN, or non-negative integers
#'
#' @param count_matrix The count matrix
#' @return TRUE if all values in count_matrix are non-NA/NaN, non-negative 
#' integers; otherwise FALSE.
areAllNonnegativeIntegers<-function(count_matrix)
{
  for (i in seq_len(nrow(count_matrix))) {
    for (j in seq_len(ncol(count_matrix))) {
      val = count_matrix[i,j]
      if(is.na(val)) {
        next # NA and NaN are allowed in probeCounts and handled later
      }
      if(!is.numeric(val)) {
        return(FALSE)
      }
      if((val < 0) | (val != round(val))) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#' Helper function to get the average volume of nucleus sampled
#' given the nucleus radius and tissue section height (thickness).
#'
#' @param r The nuclear radius
#' @param h The section height (thickness)
#' @return The average volume of nucleus sampled given the nucleus radius and 
#' section height
getAverageVolumeFrac<-function(r, h)
{
  Vavg=pi*h*(2/3*r^2 +r*h/3 - h^2/6)
  Vsphere=4/3*pi*r^3
  return(Vavg/Vsphere)
}

#' Helper function to get the maximum possible volume of nucleus
#' sampled given the nucleus radius and section height
#'
#' @param r The nuclear radius
#' @param h The section height (thickness)
#' @return The maximum possible volume of nucleus sampled given the nucleus 
#' radius and section height
getMaxVolumeFrac<-function(r,h)
{
  d=0
  hh=h/2
  v=pi*hh*(r^2 -d^2 -d*hh -(1/3)*(hh^2))
  Vmax=2*v
  Vsphere=(4/3)*pi*(r^3)
  return(Vmax/Vsphere)
}

#' Helper function that returns the minimum possible volume of nucleus
#' sampled given the nucleus radius and section height
#' 
#' @param r The radius of the nuclei
#' @param h The section height (thickness)
#' @return The minimum possible volume of nucleus sampled given the nucleus
getMinVolumeFrac<-function(r,h)
{
  d=r-h
  Vmin=pi*h*(r^2 - d^2 - d*h - (1/3)*(h^2))
  Vsphere=(4/3)*pi*(r^3)
  return(Vmin/Vsphere)
}

#' Helper function to check the arguments input to getManualCountsEstimates.
#'
#' @param probeCounts A matrix of manual spot counts with columns for 
#' probes and rows for nuclei
#' @param radius The cells' nuclear radius (must be measured in same unit as 
#' \code{height})
#' @param height The section height (must be measured in same unit as 
#' \code{radius})
#' @param volumeFracCorrection The method used to correct for volume fraction 
#' (must be 'avg', 'max', or 'min'; defaults to 'avg')
#' @return Nothing if all checks are passed; otherwise throws an error or 
#' warning message
checkManualCountsEstimatesArguments<-function(probeCounts, radius, height, 
                                              volumeFracCorrection)
{
  if(!(volumeFracCorrection %in% c("avg", "max", "min"))) {
    stop("volumeFracCorrection must be 'avg', 'max', or 'min'")
  }
  if(!is.numeric(radius)) {stop("nuclear radius must be numeric")}
  if(!is.numeric(height)) {stop("section height must be numeric")}
  if(radius <= 0) {stop("nuclear radius must be greater than 0")}
  if(height <= 0) {stop("section height must be greater than 0")}
  if(!is.matrix(probeCounts)) {stop("probeCounts must be a matrix")}
  if(ncol(probeCounts) < 1) {stop("probeCounts must have at least one column")}
  if(nrow(probeCounts) < 1) {stop("probeCounts must have at least one row")}
  if(any(is.infinite(probeCounts))) {
    stop("probeCounts cannot have any Inf or -Inf values")
  }
  if(!areAllNonnegativeIntegers(probeCounts)) {
    stop("All non-NA/NaN counts in probeCounts must be non-negative integers")
  }
  for (r in seq_len(ncol(probeCounts))) {
    if(all(is.na(probeCounts[,r]))) {
      warning("All counts of ", toString(colnames(probeCounts)[r]), 
              " probe are NA or NaN. Estimated count will be NA")
    }
    else if(all(probeCounts[,r][!is.na(probeCounts[,r])] == 0)) {
      warning("All non-NA/NaN counts of ", toString(colnames(probeCounts)[r]),
              " probe are 0. Estimated count for this ",
              "probe may not be accurate")
    }
  }
}

#' FrenchFISH function for generating volume adjusted spot counts from spots
#' which have been manually counted (uses a Markov chain Monte Carlo method).
#' Returns five quantile estimates and the mean estimate.
#'
#' @param probeCounts A matrix of manual spot counts with columns for 
#' probes and rows for nuclei
#' @param radius The cells' nuclear radius (must be measured in same unit as 
#' \code{height})
#' @param height The section height (must be measured in same unit as 
#' \code{radius})
#' @param volumeFracCorrection The method used to correct for volume fraction 
#' (must be 'avg', 'max', or 'min'; defaults to 'avg')
#' @return The volume adjusted spot counts for each probe that have been 
#' generated using MCMC modelling
#' @export
#' @examples
#' manualCountsEstimates<-getManualCountsEstimates(cbind(red=c(0,2,4),
#'     green=c(5,3,1), blue=c(3,0,2)), 8, 4)
getManualCountsEstimates<-function(probeCounts, radius, height, volumeFracCorrection = "avg")
{
  checkManualCountsEstimatesArguments(probeCounts, radius, height, volumeFracCorrection)
  
  if (volumeFracCorrection == "avg") {
    adjustFact<-getAverageVolumeFrac(radius, height)
  } else if (volumeFracCorrection == "max") {
    adjustFact<-getMaxVolumeFrac(radius, height)
  } else {
    adjustFact<-getMinVolumeFrac(radius, height)
  }

  countEstimates<-c()

  for(i in seq_len(ncol(probeCounts)))
  {
    prCts <- probeCounts[,i]
    # get probe counts column with all NA and NaN entries removed
    no_na_probeCounts <- prCts[!is.na(prCts)]

    res<-c()
    # if all counts for a probe are NA or NaN, make estimates NA
    if(length(no_na_probeCounts) == 0) {
      res<-data.frame(X2.5 = NA, X25 = NA, X50 = NA, X75 = NA, X97.5 = NA, mean = NA)
    }
    else {
      posterior_aqua<-MCMCpack::MCpoissongamma(no_na_probeCounts,
                                               0.01,0.01,5000)/adjustFact
      qtls<-summary(posterior_aqua)$quantiles
      res<-data.frame(X2.5 = c(qtls[1]), X25 = c(qtls[2]), X50 = c(qtls[3]), 
                      X75 = c(qtls[4]), X97.5 = c(qtls[5]), mean = mean(posterior_aqua))
    }
    # Append probe results to all results
    rownames(res)<-colnames(probeCounts)[i]
    countEstimates<-rbind(countEstimates, res)
  }
  return(data.frame(Probe=row.names(countEstimates), countEstimates, 
                    row.names = NULL))
}

#' Helper function to convert spot counts and nuclear area measurements
#' into continuous events for Poisson point estimation.
#'
#' @param area The nuclear area
#' @param spots The number of spots counted
#' @return Vector of continuous events for Poisson point estimation
generatePPdat<-function(area, spots)
{
  indat<-0
  for(i in seq_len(length(area)))
  {
    if(spots[i]==0)
    {
      indat[length(indat)]<-indat[length(indat)]+area[i]
    }
    else
    {
      for(k in seq_len(spots[i]))
      {
        indat<-c(indat,area[i]/spots[i])
      }
    }
  }
  return(indat)
}

#' Helper function to check the arguments input to getAutomaticCountsEstimates.
#'
#' @param probeCounts A matrix where the first column contains the areas of 
#' the nuclear blobs (this column must be named "area" and the unit of its 
#' entries must be the square of the unit used to measure \code{radius} and 
#' \code{height}) and the remaining columns (one per probe) contain the spot 
#' counts for different probes in each nuclear blob
#' @param radius The cells' nuclear radius (must be measured in same unit as 
#' \code{height})
#' @param height The section height (must be measured in same unit as 
#' \code{radius})
#' @param volumeFracCorrection The method used to correct for volume fraction 
#' (must be 'avg', 'max', or 'min'; defaults to 'avg')
#' @return Nothing if all checks are passed; otherwise throws an error or 
#' warning message
checkAutomaticCountsEstimatesArguments<-function(probeCounts, radius, height, 
                                                 volumeFracCorrection)
{
  if(!(volumeFracCorrection %in% c("avg", "max", "min"))) {
    stop("volumeFracCorrection must be 'avg', 'max', or 'min'")
  }
  if(!is.numeric(radius)) {stop("nuclear radius must be numeric")}
  if(!is.numeric(height)) {stop("section height must be numeric")}
  if(radius <= 0) {stop("nuclear radius must be greater than 0")}
  if(height <= 0) {stop("section height must be greater than 0")}
  if(!is.matrix(probeCounts)) {stop("probeCounts must be a matrix")}
  if(nrow(probeCounts) < 1) {stop("probeCounts must have at least one row")}
  if(ncol(probeCounts) < 2) {
    stop("probeCounts must have at least one column for nuclear blob ",
        "area and one column for spot counts at a probe")
  }
  if(colnames(probeCounts)[1] != "area") {
    stop('First column of probeCounts must be named "area" and contain ',
         'nuclear blob areas')
  }
  if(any(is.na(probeCounts[,1]))) {
    stop("Areas in first column cannot have any NA or NaN values")
  }
  if(any(is.infinite(probeCounts))) {
    stop("probeCounts cannot have any Inf or -Inf values")
  }
  if(!areAllNonnegativeIntegers(as.matrix(data.frame(probeCounts[,-1])))) {
    stop("All non-NA/NaN counts in probeCounts must be non-negative integers")
  }
  if(any(probeCounts[,1] <= 0)) {
    stop("All values in area column must be greater than 0")
  }
  
  for (r in 2:ncol(probeCounts)) {
    if(all(is.na(probeCounts[,r]))) {
      warning("All counts of ", toString(colnames(probeCounts)[r]), 
              " probe are NA or NaN. Estimated count will be NA")
    }
    else if(all(probeCounts[,r][!is.na(probeCounts[,r])] == 0)) {
      warning("All non-NA/NaN counts of ", toString(colnames(probeCounts)[r]),
              " probe are 0. Estimated count for this probe may not be ",
              "accurate")
    }
  }
}

#' Function to convert CSV output of the FISHalyseR automatic FISH 
#' splot counting software to a count matrix suitable for input to 
#' frenchFISH's getAutomaticCountsEstimates.
#'
#' @param pathToFishalyserCsv The path to the CSV file of automatic spot 
#' counts outputted by FISHalyseR
#' @return A count matrix suitable for input to getAutomaticCountsEstimates
#' @export
#' @examples
#' probeCounts<-convertFishalyserCsvToCountMatrix(
#'     system.file("extdata", "SampleFISH.jpg_data.csv", package="frenchFISH"))
convertFishalyserCsvToCountMatrix<-function(pathToFishalyserCsv)
{
  if (!(is.character(pathToFishalyserCsv))) {
    stop("pathToFishalyserCsv must be a string")
  }
  if (!file.exists(pathToFishalyserCsv)) {
    stop(pathToFishalyserCsv, ' does not exist')
  }

  fishalyser_df <- read.csv(pathToFishalyserCsv)
  
  # Check that fishalyser_df has the necessary columns
  if (!("area.of.nucleus" %in% colnames(fishalyser_df))) {
    stop("'area.of.nucleus' column not present in ", pathToFishalyserCsv)
  }
  
  # Get probe count columns
  count_cols <- grep("^num[.]of[.].*[.]probes$", colnames(fishalyser_df), 
                     value=TRUE)
  if (identical(count_cols, character(0))) {
    stop("No probe count columns matching regex ",
         "'^num[.]of[.].*[.]probes$' found in ", pathToFishalyserCsv)
  }
  
  # Create matrix from count columns
  fishalyser_counts_mat <- as.matrix(fishalyser_df[, count_cols])
  
  # Rename count columns to be just the color name 
  # (e.g. 'num.of.red.probes' -> 'red')
  probe_names <- unlist(strsplit(colnames(fishalyser_counts_mat), 
                                 "(num[.]of[.]|[.]probes)"))
  probe_names <- probe_names[probe_names != ""]
  colnames(fishalyser_counts_mat) <- probe_names
  
  return(cbind(area=fishalyser_df[, "area.of.nucleus"], fishalyser_counts_mat))
}

#' FrenchFISH function for generating Poisson point estimates of spot counts 
#' from spot counts which have been automatically generated. Returns a low
#' confidence interval, a median, and a high confidence interval estimate.
#'
#' @param probeCounts A matrix where the first column contains the areas of 
#' the nuclear blobs (this column must be named "area" and the unit of its 
#' entries must be the square of the unit used to measure \code{radius} and 
#' \code{height}) and the remaining columns (one per probe) contain the spot 
#' counts for different probes in each nuclear blob
#' @param radius The cells' nuclear radius (must be measured in same unit as 
#' \code{height})
#' @param height The section height (must be measured in same unit as 
#' \code{radius})
#' @param volumeFracCorrection The method used to correct for volume fraction 
#' (must be 'avg', 'max', or 'min'; defaults to 'avg')
#' @return The Poisson point estimates of spot counts for each probe
#' @export
#' @examples
#' automaticCountsEstimates<-getAutomaticCountsEstimates(
#'     cbind(area=c(250,300,450), 
#'     red=c(0,2,4), 
#'     green=c(5,3,1), 
#'     blue=c(3,0,2)), 8, 4)
getAutomaticCountsEstimates<-function(probeCounts, radius, height, 
                                      volumeFracCorrection = "avg")
{
  checkAutomaticCountsEstimatesArguments(probeCounts, radius, height,
                                         volumeFracCorrection)
  
  if (volumeFracCorrection == "avg") {
    adjustFact<-getAverageVolumeFrac(radius, height)
  } else if (volumeFracCorrection == "max") {
    adjustFact<-getMaxVolumeFrac(radius, height)
  } else {
    adjustFact<-getMinVolumeFrac(radius, height)
  }

  cellarea<-pi*(radius^2)
  countEstimates<-c()
  for(i in 2:ncol(probeCounts))
  {
    prCts <- probeCounts[,i]
    # probe counts column with all NA and NaN entries removed
    no_na_probeCounts <- prCts[!is.na(prCts)]
    # area column with all entries removed where probe counts are NA or NaN
    no_na_area <- probeCounts[,1][!is.na(prCts)]
    # if all counts for a probe are NA or NaN, make estimates NA
    if(length(no_na_probeCounts) == 0) { 
      res<-data.frame(lowCI=c(NA), median=c(NA), highCI=c(NA))
    }
    else {
      ppindat<-generatePPdat(no_na_area, no_na_probeCounts)
      # remove NA and NaA values from ppindat before cumsum
      ppcumsum<-cumsum(ppindat)/cellarea 
      PPestimate<-NHPoisson::fitPP.fun(posE=ppcumsum,nobs=max(ppcumsum),
                                       start=list(b0=0),modCI=TRUE,
                                       CIty="Transf",dplot=FALSE)
      res<-data.frame(lowCI=PPestimate@LIlambda[1],median=exp(PPestimate@coef),
                      highCI=PPestimate@UIlambda[1])
    }
    # Append probe results to all results
    rownames(res)<-colnames(probeCounts)[i]
    countEstimates<-rbind(countEstimates,res)
  }
  countEstimates<-countEstimates/adjustFact
  return(data.frame(Probe=row.names(countEstimates), countEstimates, 
                    row.names = NULL))
}
