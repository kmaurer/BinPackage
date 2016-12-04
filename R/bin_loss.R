#' Creating bin definition lists from vectors of bounds and centers
#'
#' @description Take in reduced-binned data, output loss summaries
#'
#' @param reduced_data Reduced binned data from bivariate binning function
#' @param type Type of loss summary "average" or "total"
#' @return Vector of loss summaries (both spatial and frequency losses)
#' @examples
#' red_data <- rect_bin_2d(data=diamonds, xcol="carat",ycol="price",originx=0,originy=0,widthx=.25,widthy=1000, output="reduced")
#' red_freq_data <- freq_bin(red_data,  bin_type="standard", count_type="log10", n_freq=5, pretty=FALSE, freq_col="freq")
#' bin_loss(red_data)
#' bin_loss(red_freq_data)
#' bin_loss(red_data, type="total")
bin_loss <- function(reduced_data, type="average"){
  if(type=="average"){
    loss <- with(reduced_data,
                 c(avg_loss = sum(bin_spat_loss)/sum(freq),
                   avg_standardized_loss = sum(bin_standardized_spat_loss)/sum(freq),
                   avg_loss_x = sum(bin_spat_loss_X)/sum(freq), 
                   avg_loss_y = sum(bin_spat_loss_Y)/sum(freq)))
  }
  if(type=="total"){
    loss <- with(reduced_data,
                 c(loss = sum(bin_spat_loss),
                   standardized_loss = sum(bin_standardized_spat_loss),
                   loss_x = sum(bin_spat_loss_X), 
                   loss_y = sum(bin_spat_loss_Y)))
  }
  if("freq_loss" %in% names(reduced_data)){
    loss <- c(loss,freq_loss = sum(reduced_data$freq_loss))
  }
  if("log_freq_loss" %in% names(reduced_data)){
    loss <- c(loss,log_freq_loss = sum(reduced_data$log_freq_loss))
  }
  if("log10_freq_loss" %in% names(reduced_data)){
    loss <- c(loss,log10_freq_loss = sum(reduced_data$log10_freq_loss))
  }
  return(loss)
}
