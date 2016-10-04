#------------------------------------------------------------------------------------------
### Binning functions: 
# 1d binning
  # Binning by definition (bounds and centers)
  # Creating bin definition lists from vectors of bounds and centers
  # General rectangular binning
  # Standard rectangular binning
  # Quantile rectangular binning
  # Random rectangular binning
  # Frequency binning of reduced binned data

# 2d binning
  # Independant 1d binning grid
  # Iterative conditional binning


#----------------------------------
### univariate binning based on bins definition (potentially developed from different data source, i.e. train bins then bin test)
# define bin boundaries, adapt outermost for the range of data to be binned, then use .bincode to drop to primitive to relabel
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
                       labels= labs)
  
  if(output=="centers") return(data_bins)
  if(output=="definition") return(bin_definition)
  if(output=="all") return(list(data_bins=data_bins,bin_labels=bin_labels,bin_definition=bin_definition))
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
  if(output=="all") return(list(data_bins=data_bins,bin_centers=round(bin_centers,10)))
}
# Testing with junk objects
# data.frame(xs=1:10, centers=rand_rect_bin_1d(1:10, origin=.5, width=2))


#----------------------------------
## Independent 2d Rectangular Binning for Binned Scatterplots

library(dplyr)
library(ggplot2)
xs <- diamonds$carat ; ys <- diamonds$price ; originx=0 ; originy=0 ; widthx=1 ; widthy=1000

rect_bin_2d <- function(xs, ys, originx, originy, widthx, widthy, output="full"){
  tempdat <- data.frame(xs = xs, 
                        ys=ys,
                        binxs = rect_bin_1d(xs,originx,widthx),
                        binys = rect_bin_1d(ys,originy,widthy))
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
                bin_spat_loss = sum(sqrt((xs-binxs)^2+(ys-binys)^2)),
                bin_standardized_spat_loss =  sum(sqrt((standarizedxs-standarizedbinxs)^2+(standarizedxs-standarizedbinys)^2)))
      return(reduced_binned)
  }
} 
# Test function
red_data <- rect_bin_2d(xs=diamonds$carat,ys=diamonds$price,originx=0,originy=0,widthx=1,widthy=1000, output="reduced")
red_data <- red_data %>% mutate(avg_bin_loss = bin_standardized_spat_loss/freq)
head(red_data)



#----------------------------------
## Frequency Binning  
# allows for standard or quantile binning of counts that may be raw,log,log10 (6 combinations)
# input requires binned data output, number of freq breaks and type of freq binning
# output of frequency bin values, labels and loss are attached the original then returned
freq_bin <- function(reduced_data, bin_type="standard", count_type="raw", ncolor, pretty=FALSE, freq_col="freq"){ 
  # Apply count transformations if needed
  cs <- as.data.frame(reduced_data)[,freq_col]
  if(count_type=="log")  cs <- log(cs)
  if(count_type=="log10")  cs <- log10(cs)
  # Find bin boundaries based on specification
  if(bin_type=="standard"){
    if (pretty==TRUE) {
      bounds = pretty(cs,n=ncolor)
    } else {
      round_low=floor(min(cs) / ncolor) * ncolor
      round_high=ceiling(max(cs) / ncolor) * ncolor
      bounds <- seq(round_low,round_high,by = (round_high-round_low)/ncolor )
    }
  }
  if(bin_type=="quantile"){
    quants <- quantile(cs, seq(0, 1, by=1/(2*ncolor)))
    bounds <- quants[seq(1,length(quants)+1, by=2)]
  }
  # Make bins and labels
  temp <- gen_rect_bin_1d(cs, bounds, output="all")
  reduced_data$freq_bins <- temp$data_bins
  reduced_data$freq_bin_labs <- temp$bin_labels
  
  return(reduced_data)
}
# Test out frequency binning function
# red_freq_data <- freq_bin(red_data,  bin_type="standard", count_type="raw", ncolor=5, pretty=FALSE, freq_col="freq")

head(red_freq_data)

library(RColorBrewer)
ggplot()+
  geom_tile(aes(x=binxs, y=binys, fill=as.factor(freq_bins)), data=red_freq_data) +
  scale_color_brewer()+
  theme_bw()

)))

# average standardized spatial loss = .584 standard deviations
sum(red_data$bin_standardized_spat_loss)/sum(red_data$freq)
# In context this is an average of .28 carats and $0 dollars or 
#  $500 (max loss for $)   
# in lost information for the average diamonds
sd(diamonds$carat)*.584


head(rect_bin_2d(xs=diamonds$carat,ys=diamonds$price,originx=0,originy=0,widthx=1,widthy=1000, output="centers"))
head(rect_bin_2d(xs=diamonds$carat,ys=diamonds$price,0,0,1,1000))



reduced_bin <- function(xs, ys, originx, originy, widthx, widthy, loss="spatial average")
  reduced_binned <- rect_bin_2d(xs, ys, originx, originy, widthx, widthy, output="full") %>%
    group_by(binxs,binys) %>%
    summarize(freq = n(),
              binSpatial = sum(sqrt((xs-binxs)^2+(ys-binys)^2)) )

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



