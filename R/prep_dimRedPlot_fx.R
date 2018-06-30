
#' @name fx.add_factors_to_redDim_df
#' @title Adding color values to a data.frame for plotting reduced dimensions
#' @description
#'
#' Adding color values to a data.frame for plotting reduced dimensions
#'
#' @details Basically a wrapper around \code{\link{fx.choose_vis_values_sce}}
#' that already makes smart assumptions about reasonable checks for any of
#' the factors specified via \code{size_by}, \code{color_by}, \code{shape_by}.
#' This function will add columns for those factors so that they can be
#' conveniently used for faceting etc. using \code{ggplot2}.
#'
#' @param object \code{SCE} or \code{SCESet} object
#' @param dt_to_plot \code{data.frame} fit for ggplot2-based routines
#' @param color_by entry of either \code{colData(object)},
#' \code{rowData(object)}, \code{rownames(object)},
#' which will be used to assign either discrete or continous color schemes.
#' Alternatively, you can supply a \code{list} with a label stored in \code{list$title}
#' and values stored in \code{list$result} (\code{result} should be as long as
#' \code{dim(object)[2]}. Default: \code{NULL}
#' @param shape_by default: \code{NULL}
#' @param size_by default: \code{NULL}
#' @param circle_by Used for encircling points of interest; default: \code{NULL}
#' @param exprs_values The name of the expression values which will be used for
#' the legend of the plot.
#' 
#' @return \code{data.frame} fit for plotting with \code{ggplot2} with
#' additional columns for color, shape and size
#' 
#' @seealso \code{\link{fx.choose_vis_values_sce}}
fx.add_factors_to_redDim_df <- function(object, 
    df_to_plot,
    color_by = NULL, 
    shape_by = NULL,
    size_by = NULL, 
    circle_by = NULL,
    exprs_values = NULL){
  
    df_with_cols <- df_to_plot
  
    ## Check arguments are valid and if yes, add them
    ## COLOR.........................................
    if(!is.null(color_by)){
        # user-defined label and values supplied via a list
        if(is.list(color_by)){
            ABCutilities::check_columns(c("result","title"), color_by, "list given as color_by","fx.add_factors_to_redDim_df")
            df_with_cols$color_by <- color_by$result
            names(df_with_cols)[names(df_with_cols)=="color_by"] <- color_by$title
        }else{
            color_by_out <- fx.choose_vis_values_sce(object, color_by, cell_control_default = TRUE, check_features = TRUE, exprs_values = exprs_values)
            if(!is.null(color_by_out$val)){
                df_with_cols$color_by <- color_by_out$val
                names(df_with_cols)[names(df_with_cols)=="color_by"] <- color_by_out$name # this is problematic if we want multiple instances of the coloring
            }
        }
    }
  
    ## SHAPE .........................................
    if(!is.null(shape_by) && shape_by != color_by){
        shape_by_out <- fx.choose_vis_values_sce(object, shape_by, cell_control_default = TRUE, check_features = FALSE, coerce_factor = TRUE, level_limit = 10)
        
        if(!is.null(shape_by_out$val)){
            df_with_cols$shape_by <- shape_by_out$val
            names(df_with_cols)[names(df_with_cols)=="shape_by"] <- shape_by_out$name
        }
    }
  
    ## SIZE.........................................
    if(!is.null(size_by) && !(size_by %in% c(shape_by, color_by))){
        if(is.numeric(size_by)){
            message("Numeric instances of size_by are ignored. If you want to change the
                size of all points, use the plotting function for that.")
        }else{
            size_by_out <- fx.choose_vis_values_sce(object, size_by, check_features = TRUE, check_coldata = TRUE, exprs_values = exprs_values)
            if(!is.null(size_by_out$val)){
                df_with_cols$size_by <- size_by_out$val
                names(df_with_cols)[names(df_with_cols)=="size_by"] <- size_by_out$name
            }
        }
    }
  
    ## CIRCLE.........................................
    if(!is.null(circle_by) && !(circle_by %in% c(shape_by, color_by, size_by))){
        circle_by_out <- fx.choose_vis_values_sce(object, circle_by, check_features = FALSE, check_coldata = TRUE, exprs_values = NULL, level_limit = 10)
    
        if(!is.null(circle_by_out$val)){
            df_with_cols$circle_by <- circle_by_out$val
            names(df_with_cols)[names(df_with_cols)=="circle_by"] <- circle_by_out$name
        }
    }
  
    return(df_with_cols)
}




#' @name fx.choose_vis_values_sce
#' @title Add visualization factors for ggplot2
#' @description
#'
#' Add visualization factors for ggplot2
#'
#' @details This function looks through the visualization data
#' and returns the values to be visualized. Either \code{by}
#' itself, or a column of \code{colData}, or a column of \code{rowData},
#' or the expression values of a feature.
#'
#' @param x \code{SingleCellExperiment} object
#' @param by string indicating which factor will be checked
#' @param check_coldata whether \code{colData(object)}
#' should be checked for values of \code{by}; default: TRUE
#' @param check_features whether \code{featureNames(object)} or \code{rownames(object)}
#' should be checked for values of \code{by}; default: FALSE, which saves time
#' @param exprs_values which type of expression values will be used if \code{by}
#' is a feature
#' @param level_limit maximum numbers of levels that are going to be allowed;
#' usually only applicable if \code{by} defines a shape
#'
#' @return a list with the label of the factor (\code{name}) as well as the actual
#' values that are going to be used for the plotting
#'
#' @seealso \code{\link{fx.add_factors_to_redDim_df}}
fx.choose_vis_values_sce <- function(x, by, check_coldata = TRUE,
    cell_control_default = FALSE,
    check_features = FALSE,
    exprs_values = "counts",
    coerce_factor = FALSE, 
    level_limit = NA){
  
    ## BY == character .........................................
    if (is.character(by)) {
        if (length(by) != 1L) {
            warning("'by' should be a character vector of length 1")
            vals <- NULL
        }
    
        ## checking if it refers to a field in colData/rowData.
        vals <- NULL
        if (check_coldata) {
            if (by %in% names(colData(x)) ) {
                vals <- colData(x)[[by]]
            }
            unfound <- "names(colData(x))"
        } else {
            if (by %in% names(rowData(x))) {
                vals <- rowData(x)[[by]]
            }
            unfound <- "names(rowData(x))"
        }
      
        ## checking if it refers to a feature
        if (check_features) {
            if (is.null(vals) && by %in% rownames(x)) {
                if(is.null(exprs_values)){
                    warning("You did not specify which type of expression values you wanted
                        to extract for the feature. Counts will be used per default.")
                }
                vals <- assay(x[by,], exprs_values)[,]
            }
            unfound <- append(unfound, "rownames(x)")
        }  
    
        ## throwing an error if still unfound
        if (is.null(vals)) {
            warning(sprintf("'%s' not found %s", by, paste(sprintf("in '%s'", unfound), collapse = " or ")))
            vals <- NULL
        }
    
    } else if (is.data.frame(by)) { ## BY == data.frame...........................
        if (ncol(by) != 1L) {
            stop("'by' should be a data frame with one column")
        } else if (nrow(by) != ncol(x)) {
            stop("'nrow(by)' should be equal to number of columns in 'x'")
        }
      
        ## Allow arbitrary values to be specified.
        vals <- by[,1]
        by <- colnames(by)
      
    } else { ## BY != character | BY != data.frame ==> ERROR....................
        if (!is.null(by)){
            warning("Invalid value of 'by' supplied (neither character nor data.frame)")
            vals <- NULL
        }
      
        ## Switching to cell controls if desired.................................
        if (cell_control_default) {
            if ( "is_cell_control" %in% names(colData(x)) ) {
                by <- "is_cell_control"
                vals <- colData(x)[[by]]
            }
        }
    }
  
    ## Checking the level limit..........................................
    if (coerce_factor && !is.null(vals)) {
        vals <- factor(vals)
        if (level_limit < nlevels(vals)){
            warning(sprintf("The number of unique levels exceeds %i", level_limit))
            vals <- NULL
        }
    }
  
    return(list(name = by, val = vals))

}




#' @name fx.add_cellInfo_to_redDim_df
#' @title Add metadata about the cells (columns) to the data.frame for plotting
#' @description
#'
#' Add metadata about the cells (columns) to the data.frame for plotting
#'
#' @param object \code{SingleCellExperiment} object
#' @param df_to_plot \code{data.frame} for ggplot2-based plotting
#' @param what_info vector of names that will be used to extract the
#' corresponding values from \code{object}
#' @param max_levels indicate a number if you want to limit the number of levels
#' that any of the entries specified by \code{what_info} can have. If the
#' number is exceeded, the entry will be dropped. Default: NULL
#' @return \code{df_to_plot} with additional columns
#' @examples \dontrun{#' 
#' df_to_plot <- fx.add_cellInfo_to_redDim_df(sce, df_to_plot, names(colData(sce))[1:4])
#'}
#'
fx.add_cellInfo_to_redDim_df <- function(object, df_to_plot, what_info, max_levels = NULL){
  
    ## avoid duplication
    what_info <- what_info[ !what_info %in% names(df_to_plot) ]
  
    if(length(what_info) == 0){
        return(df_to_plot)
    }else{
        ## avoid error due to missing entry
        if( !all(what_info %in% names(colData(object))) ){
            missing_entry <-  what_info[ !what_info %in% names(colData(object)) ]
            warning(paste( paste(missing_entry, collapse = ","), "is/are not in colData(object). Will be dropped."))
            what_info <- what_info[ what_info != missing_entry ]
        }
  
        ## extract values
        info_df <- colData(object)[what_info]
  
        ## if wanted, check the number of levels
        if(!is.null(max_levels)){
            keep <- unlist(lapply(info_df, function(x) nlevels(as.factor(x)) <= max_levels))
            info_df <- info_df[keep]
        }
  
        ## combine
        df_with_info <- cbind(df_to_plot, data.frame(info_df))
  
        return(df_with_info)
    }
}




#' @name fx.return_axis_length
#' @title Get the number of coordinates stored in the axis value vector
#' @description
#'
#' Get the number of coordinates stored in the axis value vector
#' 
#' @description The axis coordinates for tSNE and PCA can be given either as a
#' one-dimensional vector, or as a matrix or data.frame with one column. This
#' function returns the correct length that can be used to test against the 
#' number of cells stored in the singleCellExperiment object.
#' 
#' @param whatever is going to be passed on to get_reducedDimPlot.sce as the X
#' and Y coordinate vectors
#' 
#' @return a number
#' 
#' @examples
#' 
#' df <- data.frame(a = c(1:10))
#' fx.return_axis_length(df[,1,drop=FALSE])
#' fx.return_axis_length(df[,1])
#' 
fx.return_axis_length = function(axs){
  
    if( is.null(dim(axs)) ){
        nout <- length(axs)
    }else{
        nout <- dim(axs)[1]
    } 

    return(nout)

}




#' @name fx.get_percent_var 
#' @title Extract the percent variation values stored as attributes
#' @description
#'
#' Extract the percent variation values stored as attributes
#' 
#' @param rd PCA result object
#' @param pcs_select the number of PCs for which the values should be returned
fx.get_percent_var = function (rd, pcs_select = NULL) {
    var_vals <- attr(rd, "percentVar")
    if (is.null(var_vals)) {
        return(var_vals)
    }else {
        names(var_vals) <- colnames(rd)
        if (!is.null(pcs_select)) {
            var_vals <- var_vals[pcs_select]
        }
        return(var_vals)
    }
}
