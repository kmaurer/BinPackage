#------------------------------------------------------------------------------------------
### Binning functions: 
# 1d binning
  # Binning by definition (bounds and centers)
  # Creating bin definition lists from vectors of bounds and centers
  # General rectangular binning
  # Standard rectangular binning
  # Quantile rectangular binning
  # Random rectangular binning

# 2d binning
  # Independant 1d binning grid
  # Iterative conditional binning


#----------------------------------
### univariate binning based on bins definition (potentially developed from different data source, i.e. train bins then bin test)
# define bin boundaries, adapt outermost for the range of data to be binned, then use .bincode to drop to primitive to relabel
bin_by_def <- function(new_data_vec, bin_definition, output="centers"){
  bounds <- bin_definition$bin_bounds
  bounds[1] <- min(min(new_data_vec),bounds[1])
  bounds[length(bounds)] <- max(max(new_data_vec), bounds[length(bounds)])
  data_bins <- bin_definition$bin_centers[.bincode(new_data_vec,bounds,right=FALSE, include.lowest=TRUE)]
  names(data_bins) <- names(new_data_vec)
  if(output=="centers") return(data_bins)
  if(output=="definition") return(bin_definition)
  if(output=="both") return(list(data_bins=data_bins,bin_definition=bin_definition))
}
# Testing with junk objects
# bin_definition<-  make_bin_def(bin_centers,bin_bounds)
# new_data_vec <- -1:11
# cbind(new_data_vec, bin_by_def(new_data_vec, bin_definition))

#----------------------------------
## Making bin definitions
make_bin_def <- function(centers,bounds) return(list(bin_centers=centers,bin_bounds=bounds))
#make_bin_def(1:10,seq(.5,11.5,by=1))

#----------------------------------
## Generalized Rectangular 1d Binning : given the set of bin breaks
gen_rect_bin_1d <- function(xs, bounds, output="centers"){
  bin_bounds <- bounds
  bin_centers <-bounds[1:(length(bounds)-1)] + (bounds[2:length(bounds)]-bounds[1:(length(bounds)-1)])/2
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}
# Testing with junk objects
# xs <- runif(100,1,25)
# bounds <- c(1,4,7,12,13,19,23,24,25,29)
# head(data.frame(xs, centers=gen_rect_bin_1d(xs, bounds ,output="centers")))

#----------------------------------
## Standard Rectangular 1d Binning: given origin and consistent bin width
rect_bin_1d <- function(xs, origin, width, output="centers"){
  bin_bounds <- origin + width*(0:(ceiling(diff(range(xs))/width)))
  bin_centers <- origin + width*(1:( ceiling(diff(range(xs))/width)) - 0.5)
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}
# Testing with junk objects
# xs <- ggplot2::diamonds$price; nbins=4; origin=min(xs); width=diff(range(xs))/nbins
# rect_bin_1d(rnorm(1000,0,1),origin=-4,width=1,output="both")

#----------------------------------
## Quantile 1d Binning
# used for binning the counts values by quantile
# define vector of counts and number of bin
quant_bin_1d <- function(xs, nbin, output="centers"){
  quants <- quantile(xs, seq(0, 1, by=1/(2*nbin)))
  bin_centers <- quants[seq(2,length(quants)-1, by=2)]
  bin_bounds <- quants[seq(1,length(quants)+1, by=2)]
  return(bin_by_def(xs, make_bin_def(bin_centers,bin_bounds),output=output))
}
# Testing with junk objects
# quant_bin_1d(ggplot2::diamonds$price,4,output="both")

#----------------------------------
## Random 1d Rectangular Binning
# Establish bounding bin centers for each observation
# Then use Unif(0,1) draws compared to assignment probs to allocate
# Reassignment for values below first center or above highest center
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
  if(output=="centers") return(data_bins)
  if(output=="definition") return(list(bin_centers=round(bin_centers,10)))
  if(output=="both") return(list(data_bins=data_bins,bin_centers=round(bin_centers,10)))
}
# Testing with junk objects
# data.frame(xs=1:10, centers=rand_rect_bin_1d(1:10, origin=.5, width=2))

#----------------------------------
## Independent 2d Rectangular Binning for Binned Scatterplots

library(dplyr)
library(ggplot2)
xs <- diamonds$carat ; ys <- diamonds$price ; originx=0 ; originy=0 ; widthx=1 ; widthy=1000


rect_bin_2d <- function(xs,ys, originx, originy, widthx, widthy, type="standard"){
  tempdat <- data.frame(xs = xs, 
                        ys=ys,
                        binxs = rect_bin_1d(xs,originx,widthx),
                        binys = rect_bin_1d(ys,originy,widthy)) 
  
  reduced_binned <- tempdat %>%
    group_by(binxs,binys) %>%
    summarize(binfreq = n(),
              binspatialloss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)) )

  summarydata <- data.frame( originx = originx, originy = originy, 
                             widthx = widthx, widthy = widthy,
                             totalSpatialLoss = sum(outdat$binspatialloss))
  templist <- list(bindat = outdat[,1:3],
                   summarydata)
  return(templist)
}


rect_bin_2d <- function(xs,ys, originx, originy, widthx, widthy, type="standard"){
  tempdat <- data.frame(xs = xs, ys=ys,
                        binxs = StandRectBin1d(xs,originx,widthx),
                        binys = StandRectBin1d(ys,originy,widthy))
  
  if(type=="random"){
    tempdat<- data.frame( tempdat,
                          randbinxs = RandRectBin1d(xs,originx,widthx),
                          randbinys = RandRectBin1d(ys,originy,widthy))
    tempdat$binname <- paste(tempdat$binxs,tempdat$binys)
    tempdat$randbinname <- paste(tempdat$randbinxs,tempdat$randbinys)
    tempdat$randdist <- with(tempdat,sqrt((xs-randbinxs)^2 + (ys-randbinys)^2))
    tempdat$index <- 1:length(xs)
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
  }
  outdat <- ddply(tempdat, .(binxs,binys), summarise,
                  binfreq = length(xs),
                  binspatialloss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)) )
  summarydata <- data.frame( originx = originx, originy = originy, 
                             widthx = widthx, widthy = widthy,
                             totalSpatialLoss = sum(outdat$binspatialloss))
  templist <- list(bindat = outdat[,1:3],
                   summarydata)
  return(templist)
}



