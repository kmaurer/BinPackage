#' Creating bin definition lists from vectors of bounds and centers
#'
#' @description Construct binned scatterplot from reduced binned data, build upon \code{geom_tile} from \code{ggplot2} 
#' Takes in either reduced-binned data with raw counts or with frequency-binned counts. Allow for raw/log/log10 scaling
#'
#' @param reduced_data Reduced binned data from either \code{rect_bin_2d} or \code{rand_rect_bin_2d}
#' @param count_scale Transformation for shading scale for bin counts
#' @param freq_binned TRUE/FALSE does reduced binned data has been frequency binned.
#' @return Vector of loss summaries
#' @examples
#' n_freq <-5
#' red_data <- rect_bin_2d(data=diamonds, xcol="carat",ycol="price",originx=0,originy=0,widthx=.25,widthy=1000, output="reduced")
#' red_freq_data <- freq_bin(red_data,  bin_type="standard", count_type="log10", n_freq=n_freq, pretty=FALSE, freq_col="freq")
#' binned_scatter(red_freq_data, count_scale="log10", freq_binned = FALSE)
#' binned_scatter(red_freq_data, count_scale="log10",freq_binned = TRUE) +
#'  xlab("Carat Weight") + ylab("Price (dollars)") +
#'  ggtitle("Binned Scatterplot of Diamond Price By Carat Weight")
binned_scatter <- function(reduced_data, count_scale="identity", freq_binned=FALSE){
  color_lab <- "Frequency"
  if(count_scale == "log") color_lab <- "Frequency (Log-Scaled)"
  if(count_scale == "log10") color_lab <- "Frequency (Log10-Scaled)"
  
  if(freq_binned==FALSE){
    p1 <- ggplot()+
      geom_tile(aes(x=binxs, y=binys, fill=freq), data=reduced_data) +
      scale_fill_gradient(color_lab,low="#56B1F7", high="#132B43", trans=count_scale) +
      theme_bw() +
      theme(legend.position="bottom", 
            legend.key.width = unit(3, "cm")) 
  }
  
  if(freq_binned==TRUE){
    nf <- length(levels(reduced_data$freq_bins))
    p1 <- ggplot()+
      geom_tile(aes(x=binxs, y=binys, fill=freq_bins), data=reduced_data) +
      scale_fill_manual(paste("Binned",color_lab),
                        guide = guide_legend(nrow=1,
                                             keywidth=1.2,
                                             keyheight=0.5,
                                             default.unit="inch",
                                             label.position="bottom",
                                             label.hjust=.5,  
                                             title.position = "top"),
                        values=seq_gradient_pal("#56B1F7", "#132B43")((1:nf)/nf),
                        labels=levels(reduced_data$freq_bin_labs),
                        drop=FALSE) +
      theme_bw() +
      theme(legend.key.width = unit(nf*.8, "cm"),
            legend.position="bottom") 
  }
  return(p1)
}


#' Creating bin definition lists from vectors of bounds and centers
#'
#' @description All-in-one Binning and Scatterplot Creating Super Function
#' Need to specify all origins, widths, transformations, frequency bins
#'
#' @param data data frame with at least two numeric columns to be binned
#' @param xcol name of column to bin for x dimension
#' @param ycol name of column to bin for y dimension
#' @param originx origin numeric value at which to begin standard rectangular binning for x dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthx Positive numeric value at which to space standard rectangular bins for x dimension.
#' @param originy origin numeric value at which to begin standard rectangular binning for y dimension. This should be selected at or below minimum value in \code{xs}.
#' @param widthy Positive numeric value at which to space standard rectangular bins for y dimension.
#' @param freq_binned TRUE/FALSE for using frequency binning for discretized color scale.
#' @param freq_bin_type "standard" or "quantile" binning specification
#' @param count_type "raw","log", or "log10" transformation of counts specified
#' @param n_freq positive integer specifying number of bins
#' @param pretty Logical for using pretty function for cleaner frequency bin breaks
#' @return Vector of loss summaries
#' @examples
#' binned_scatter_raw(data=diamonds, xcol="carat", ycol="price", originx=0,
#' originy=0, widthx=.25, widthy=1000, count_type = "log10")
#' binned_scatter_raw(data=diamonds, xcol="carat", ycol="price", originx=0,
#'                    originy=0, widthx=.25, widthy=1000, count_type="raw",
#'                    freq_binned=TRUE, freq_bin_type="quantile", n_freq=5, pretty=FALSE)
binned_scatter_raw <- function(data, xcol, ycol, originx, originy, widthx, widthy, count_type="raw",
                               freq_binned=FALSE, freq_bin_type="standard",  n_freq=1, pretty=FALSE){
  red_data <- rect_bin_2d(data=data, xcol=xcol, ycol=ycol, originx=originx,originy=originy,widthx=widthx,widthy=widthy, output="reduced")
  red_freq_data <- freq_bin(red_data,  bin_type=freq_bin_type, count_type=count_type, n_freq=n_freq, pretty=pretty, freq_col="freq")
  count_scale <- "identity" 
  if(count_type != "raw") count_scale <- count_type
  return(binned_scatter(red_freq_data, count_scale=count_scale, freq_binned = freq_binned))
}

