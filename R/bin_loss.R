#' Creating bin definition lists from vectors of bounds and centers
#'
#' @description Take in reduced-binned data, output loss summaries
#'
#' @param reduced_data Reduced binned data from bivariate binning function
#' @param type Type of loss summary "average" or "total"
#' @return Vector of loss summaries
#' @examples
#' reduced_data <- rect_bin_2d(data=diamonds, xcol="carat", ycol="price",originx=0,originy=0,widthx=.25,widthy=100, output="reduced")
#' bin_loss(reduced_data)
#' bin_loss(reduced_data, type="total")
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
  return(loss)
}
