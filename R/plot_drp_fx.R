#' Extract factor information
#' 
#' @param drp_object DimRedPlot object, ideally generated with \code{\link{}}
#' @param fctr the value of the given factor, e.g. "ENSMUSG00000025920" or 
#' "condition"
#' @param fctr_type the type of factor, i.e. the name such as "color_by", or
#' "size_by"
#' @param ignore_drp_labels either TRUE/FALSE or a vector of fctr_type(s) to be
#' ignored from the drp object
#' 
#' @return the string that will be used to assign the fct
fx.set_factors <- function(drp_object, fctr = color_by, fctr_type = "color_by",
                           ignore_drp_labels){
 
   fct_out <- fctr
  
  if(is.null(fctr)){ # else: fct_out <- fctr
    if( all(ignore_drp_labels == FALSE) | !(fctr_type %in% ignore_drp_labels) ){
      fct_out <- drp_object$factors[[fctr_type]]
    }
  }
  
  return(fct_out)
}
