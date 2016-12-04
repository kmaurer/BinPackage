#' Creating bin definition lists from vectors of bounds and centers
#'
#' @description Helper function to combine vectors of bin centers and bounds into a list to be used in other functions that bin by definition.
#'
#' @param centers vector of bin centers
#' @param bounds vector of bin boundaries
#' @return Returns list of bin definition
#' @examples
#' make_bin_def(1:10,seq(.5,11.5,by=1))
make_bin_def <- function(centers,bounds){
  return(list(bin_centers=centers,bin_bounds=bounds))
}

#' Binning by definition (bounds and centers)
#'
#' @description univariate binning based on bins definition (potentially developed from different data source, i.e. train bins then bin test)
#' define bin boundaries, adapt outermost for the range of data to be binned, then use .bincode to drop to primitive to relabel
#'
#' @param new_data_vec data vector to bin using defined partion.
#' @param bin_definition list containing bounds of binning partition.
#' @param output "centers" or "all" to select if function will return only bin centers vector or list with bin centers and defining bounds, respectively.
#' @return returns either the new binned data, or a list of the data and defining bounds
bin_by_def <- function(new_data_vec, bin_definition, output="centers"){
  # make bounds
  bounds <- bin_definition$bin_bounds
  bounds[1] <- min(min(new_data_vec),bounds[1])
  bounds[length(bounds)] <- max(max(new_data_vec), bounds[length(bounds)])
  # allocated to centers
  data_bins <- bin_definition$bin_centers[.bincode(new_data_vec,bounds,right=FALSE, include.lowest=TRUE)]
  names(data_bins) <- names(new_data_vec)
  # Add bin labels
  labs <-  paste(c("[",rep("(",(length(bin_definition$bin_bounds)-2))),
                 bin_definition$bin_bounds[1:(length(bin_definition$bin_bounds)-1)],
                 ",",bin_definition$bin_bounds[2:length(bin_definition$bin_bounds)],
                 "]",sep="")
  bin_labels <- factor(.bincode(new_data_vec,bounds,right=FALSE,include.lowest=TRUE),
                       levels= 1:(length(bounds)-1),
                       labels= labs)
  if(output=="centers") return(data_bins)
  if(output=="all") return(list(data_bins=data_bins,bin_labels=bin_labels,bin_definition=bin_definition))
}


#' General rectangular binning
#'
#' @description Generalized Rectangular 1d Binning : given the set of bin breaks
#'
#' @param xs Vector of numeric values to bin
#' @param bounds Vector of unique, increasing boundaries values that specify bin partitions
#' @param output Parameter reserved for expanding output types in future package versions.
#' @return bin centers as a vector corresponding to \code{xs} vector
#' @examples
#' xs <- runif(100,1,25)
#' bounds <- c(1,4,7,12,13,19,23,24,25,29)
#' gen_rect_bin_1d(xs, bounds ,output="centers")
gen_rect_bin_1d <- function(xs, bounds, output="centers"){
  bin_bounds <- bounds
  bin_centers <-bounds[1:(length(bounds)-1)] + (bounds[2:length(bounds)]-bounds[1:(length(bounds)-1)])/2
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}

#' Standard rectangular binning
#'
#' @description Standard Rectangular 1d Binning: given origin and consistent bin width.
#'
#' @param xs Vector of numeric values to bin.
#' @param origin numeric value at which to begin standard rectangular binning. This should be selected at or below minimum value in \code{xs}.
#' @param width Positive numeric value at which to space standard rectangular bins.
#' @param output Parameter reserved for expanding output types in future package versions.
#' @return bin centers as a vector corresponding to \code{xs} vector
#' @examples
#' xs <- ggplot2::diamonds$price; nbins=4; origin=min(xs); width=diff(range(xs))/nbins
#' rect_bin_1d(rnorm(1000,0,1),origin=-4,width=1)
rect_bin_1d <- function(xs, origin, width, output="centers"){
  bin_bounds <- origin + width*(0:(ceiling(diff(range(xs))/width)))
  bin_centers <- origin + width*(1:( ceiling(diff(range(xs))/width)) - 0.5)
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}


#' Quantile rectangular binning
#'
#' @description Quantile 1d Binning: used for binning the counts values by quantile, must define vector of counts and number of bins.
#'
#' @param xs Vector of numeric values to bin.
#' @param nbin Number of quantile based bins, must be a positive integer.
#' @param output Parameter reserved for expanding output types in future package versions.
#' @return bin centers as a vector corresponding to \code{xs} vector
#' @examples
#' quant_bin_1d(ggplot2::diamonds$price,4)
quant_bin_1d <- function(xs, nbin, output="centers"){
  quants <- quantile(xs, seq(0, 1, by=1/(2*nbin)))
  bin_centers <- quants[seq(2,length(quants)-1, by=2)]
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}


#' Random 1d Rectangular Binning
#'
#' @description Establish bounding bin centers for each observation
#' Then use Unif(0,1) draws compared to assignment probs to allocate
#' Reassignment for values below first center or above highest center
#'
#' @param xs Vector of numeric values to bin.
#' @param origin numeric value at which to begin standard rectangular binning. This should be selected at or below minimum value in \code{xs}.
#' @param width Positive numeric value at which to space standard rectangular bins.
#' @param output Parameter reserved for expanding output types in future package versions.
#' @return bin centers as a vector corresponding to \code{xs} vector
#' @examples
#' rand_rect_bin_1d(1:10, origin=.5, width=2)
rand_rect_bin_1d <- function(xs, origin, width, output="centers"){
  bin_centers <- origin + width*(1:( ceiling(diff(range(xs))/width)) - 0.5)
  #create sets of lower and upper bounding centers
  lbs <- bin_centers[1] + width*floor((xs-bin_centers[1])/width)
  ubs <- bin_centers[1] + width*ceiling((xs-bin_centers[1])/width)
  # initially assign all values to upper bound
  data_bins <- ubs
  # then use runif to reassign based on distance from
  plower <- (ubs - xs)/width
  lowerindex <- (plower > runif(length(xs), 0, 1))
  data_bins[lowerindex] <- lbs[lowerindex]
  data_bins[xs < bin_centers[1]] <- bin_centers[1]
  data_bins[xs > bin_centers[length(bin_centers)]] <- bin_centers[length(bin_centers)]
  # Return output based on option selected
  return(data_bins)
}


#' Frequency binning of reduced binned data
#'
#' @description allows for standard or quantile binning of bin counts that may be raw,log,log10 (6 combinations)
#' for use in discretizing binned scatterplot shading.
#'
#' @param reduced_data Requires reduced binned data output from \code{rect_bin_2d} function, number of freq breaks and type of freq binning
#' @param bin_type "standard" or "quantile" binning specification
#' @param count_type "raw","log", or "log10" transformation of counts specified
#' @param n_freq positive integer specifying number of bins
#' @param pretty Logical for using pretty function for cleaner frequency bin breaks
#' @param freq_col name of frequency column to be binned
#' @return output of frequency bin values, labels and loss are attached the original then returned
#' @examples
#' reduced_data <- rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="reduced")
#' freq_red_data <- freq_bin(reduced_data, bin_type="standard", count_type="raw", n_freq=5, pretty=FALSE, freq_col="freq")
#' head(freq_red_data)
#' freq_log10_red_data <- freq_bin(reduced_data, bin_type="standard", count_type="log10", n_freq=5, pretty=FALSE, freq_col="freq")
#' head(freq_log10_red_data)
#' sum(freq_log10_red_data$log10_freq_loss)
freq_bin <- function(reduced_data, bin_type="standard", count_type="raw", n_freq=5, pretty=FALSE, freq_col="freq"){
  # Apply count transformations if needed
  cs <- as.data.frame(reduced_data)[,freq_col]
  if(count_type=="log")  cs <- log(cs)
  if(count_type=="log10")  cs <- log10(cs)
  # Find bin boundaries based on specification
  if(bin_type=="standard"){
    if (pretty==TRUE) {
      bounds = pretty(cs,n=n_freq)
    } else {
      round_low=floor(min(cs) / n_freq) * n_freq
      round_high=ceiling(max(cs) / n_freq) * n_freq
      bounds <- seq(round_low,round_high,by = (round_high-round_low)/n_freq )
    }
  }
  if(bin_type=="quantile"){
    quants <- quantile(cs, seq(0, 1, by=1/(2*n_freq)))
    bounds <- quants[seq(1,length(quants)+1, by=2)]
  }
  # Make bins and labels
  temp <- gen_rect_bin_1d(cs, bounds, output="all")
  if(count_type=="raw") reduced_data$freq_loss <- (cs-temp$data_bins)^2
  if(count_type=="log"){
    reduced_data$log_freq <- cs
    reduced_data$log_freq_loss <- (cs-temp$data_bins)^2
  }
  if(count_type=="log10"){
    reduced_data$log10_freq <- cs
    reduced_data$log10_freq_loss <- (cs-temp$data_bins)^2
  }
  reduced_data$freq_bins <- factor(temp$data_bins, levels=temp$bin_definition$bin_centers)
  reduced_data$freq_bin_labs <- temp$bin_labels
  return(reduced_data)
}


#' Independant 2d binning grid
#'
#' @description Independent 2d Rectangular Binning for Binned Scatterplots. Built to return a data frame containing:
#' full data - the input data with binned xs and ys as new columns,
#' OR
#' centers data - the new binned xs and ys centers only,
#' OR
#' reduced binned data - one row per bivariate bin with columns for x and y centers, bin counts, and euclidean spatial losses per bin.
#'
#' @param data data frame with at least two numeric columns to be binned
#' @param xcol name of column to bin for x dimension
#' @param ycol name of column to bin for y dimension
#' @param originx origin numeric value at which to begin standard rectangular binning for x dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthx Positive numeric value at which to space standard rectangular bins for x dimension.
#' @param originy origin numeric value at which to begin standard rectangular binning for y dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthy Positive numeric value at which to space standard rectangular bins for y dimension.
#' @param output "full" (default),"centers",or"reduced"
#' @return returns data frame of specified type
#' @examples
#' full_data <- rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="full")
#' head(full_data)
#' centers_data <- rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="centers")
#' head(centers_data)
#' reduced_data <- rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="reduced")
#' head(reduced_data)
rect_bin_2d <- function(data, xcol, ycol, originx, originy, widthx, widthy, output="full"){
  tempdat <- data.frame(xs = data.frame(data)[,xcol],
                        ys=data.frame(data)[,ycol],
                        binxs = rect_bin_1d(data.frame(data)[,xcol],originx,widthx),
                        binys = rect_bin_1d(data.frame(data)[,ycol],originy,widthy))
  if(output=="full") return(tempdat)
  if(output=="centers") return(tempdat[,c("binxs","binys")])
  if(output=="reduced") {
    reduced_binned <- tempdat %>%
      mutate(standarizedxs = (xs-mean(xs))/sd(xs),
             standarizedys = (ys-mean(ys))/sd(ys),
             standarizedbinxs = (binxs-mean(xs))/sd(xs),
             standarizedbinys = (binys-mean(ys))/sd(ys)) %>%
      group_by(binxs,binys) %>%
      summarize(freq = n(),
                bin_spat_loss_X = sum(abs(xs-binxs)),
                bin_spat_loss_Y = sum(abs(ys-binys)),
                bin_spat_loss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)),
                bin_standardized_spat_loss =  sum(sqrt((standarizedxs-standarizedbinxs)^2+(standarizedxs-standarizedbinys)^2)))
    return(data.frame(reduced_binned))
  }
}

#' 2d Random Rectangular Binning for Binned Scatterplots
#'
#' @description 2d Random Rectangular Binning for Binned Scatterplots. Built to return a data frame containing:
#' full data - the input data with binned xs and ys as new columns,
#' OR
#' centers data - the new binned xs and ys centers only,
#' OR
#' reduced binned data - one row per bivariate bin with columns for x and y centers, bin counts, and net euclidean spatial losses per bin (after post-processing to correct for perceived loss).
#'
#' @param data data frame with at least two numeric columns to be binned
#' @param xcol name of column to bin for x dimension
#' @param ycol name of column to bin for y dimension
#' @param originx origin numeric value at which to begin random rectangular binning for x dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthx Positive numeric value at which to space random rectangular bins for x dimension.
#' @param originy origin numeric value at which to begin random rectangular binning for y dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthy Positive numeric value at which to space random rectangular bins for y dimension.
#' @param output "full" (default),"centers",or"reduced"
#' @return returns data frame of specified type
#' @examples
#' full_data <- rand_rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="full")
#' head(full_data)
#' centers_data <- rand_rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="centers")
#' head(centers_data)
#' reduced_data <- rand_rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.5,widthy=1000, output="reduced")
#' head(reduced_data)
rand_rect_bin_2d <- function(data, xcol, ycol, originx, originy, widthx, widthy, output="full"){
  tempdat <- data.frame(xs = data.frame(data)[,xcol], 
                        ys=data.frame(data)[,ycol],
                        randbinxs = rand_rect_bin_1d(data.frame(data)[,xcol],originx,widthx),
                        randbinys = rand_rect_bin_1d(data.frame(data)[,ycol],originy,widthy),
                        binxs = rect_bin_1d(data.frame(data)[,xcol],originx,widthx),
                        binys = rect_bin_1d(data.frame(data)[,ycol],originy,widthy))
  tempdat$binname <- paste(tempdat$binxs,tempdat$binys)
  tempdat$randbinname <- paste(tempdat$randbinxs,tempdat$randbinys)
  tempdat$randdist <- with(tempdat,sqrt((xs-randbinxs)^2 + (ys-randbinys)^2))
  tempdat$index <- 1:nrow(tempdat)
  # points where mismatch between standard and random binning
  mismatchindex <- which(tempdat$binxs != tempdat$randbinxs | tempdat$binys != tempdat$randbinys )
  mmdat <- tempdat[mismatchindex,]
  # identify which need to be swapped to standard for net spatial loss
  #loop over all neighboring bin pairs (shift all xs over by 1 then all )
  xbins <- seq(min(c(tempdat$binxs,tempdat$randbinxs)) , max(c(tempdat$binxs,tempdat$randbinxs)), by=widthx)
  ybins <- seq(min(c(tempdat$binys,tempdat$randbinys)) , max(c(tempdat$binys,tempdat$randbinys)), by=widthy)
  nbrs <- data.frame(binxs = rep(rep(xbins,length(ybins)),2) ,
                     binys = rep(rep(ybins,each=length(xbins)),2),
                     nbrsxs = c(rep(xbins+widthx,length(ybins)),rep(xbins,length(ybins))),
                     nbrsys = c(rep(ybins,each=length(xbins)), rep(ybins+widthy,each=length(xbins))) )
  nbrs$binname <- paste(nbrs$binxs,nbrs$binys)
  nbrs$nbrsname <- paste(nbrs$nbrsxs,nbrs$nbrsys)   
  swapindex <- NULL    
  for (i in 1:nrow(nbrs)){
    #id points in standard bin i assigned to bin j
    itoj <- which(mmdat$binname == nbrs$binname[i] & mmdat$randbinname == nbrs$nbrsname[i])
    #id points in standard bin j assigned to bin i
    jtoi <- which(mmdat$binname == nbrs$nbrsname[i] & mmdat$randbinname == nbrs$binname[i])
    # number to swap in bins i and j is equal to minimum misplaced
    nswap <- min(length(itoj), length(jtoi))
    # if there are points to swap, then pick the ones with largest distance
    # from point to random bin center
    if(nswap > 0){ 
      swapindex <- c(swapindex,mmdat$index[itoj][order(mmdat$randdist[itoj],decreasing=TRUE)[nswap]])
      swapindex <- c(swapindex,mmdat$index[jtoi][order(mmdat$randdist[jtoi],decreasing=TRUE)[nswap]])
    }
  }
  swapindex <- unique(swapindex)
  tempdat$binxs[!(tempdat$index %in% swapindex)] <- tempdat$randbinxs[!(tempdat$index %in% swapindex)]
  tempdat$binys[!(tempdat$index %in% swapindex)] <- tempdat$randbinys[!(tempdat$index %in% swapindex)]
  if(output=="full") return(tempdat)
  if(output=="centers") return(tempdat[,c("binxs","binys")])
  if(output=="reduced") {
    reduced_binned <- tempdat %>%
      mutate(standarizedxs = (xs-mean(xs))/sd(xs),
             standarizedys = (ys-mean(ys))/sd(ys),
             standarizedbinxs = (binxs-mean(xs))/sd(xs),
             standarizedbinys = (binys-mean(ys))/sd(ys)) %>%
      group_by(binxs,binys) %>%
      summarize(freq = n(),
                bin_spat_loss_X = sum(abs(xs-binxs)),
                bin_spat_loss_Y = sum(abs(ys-binys)),
                bin_spat_loss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)),
                bin_standardized_spat_loss =  sum(sqrt((standarizedxs-standarizedbinxs)^2+(standarizedxs-standarizedbinys)^2)))
    return(data.frame(reduced_binned))
  }
} 

