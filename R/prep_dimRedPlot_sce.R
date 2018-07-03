
#' @name get_reducedDimPlot.sce
#' @title Prepare the "DimRedPlot" object for plotting reduced dimension information
#' @description
#' 
#' Wrapper function for extracting the X and Y coordinates for the
#' plot as well as the variables for coloring, shape assignment etc.
#' This function expects a SingleCellExperiment object.
#' 
#' @details
#' \emph{which_reddim} The reducedDims slot of the SingleCellExperiment object
#' is a list (\code{\link{SimpleList}}. To know which coordinates to extract,
#' the accessor for the specific reducedDims list member must be given.
#'
#' @param object \code{SingleCellExperiment} object
#' @param which_reddim accessor to retrieve the list of choice for the reduced
#' dimension coordinates from \code{reducedDims(object)}. See \code{details}.
#' @param which_pcs vector of integers to specify which (principal) components 
#' are going to be used for the x and y coordinates. The default (\code{c(1,2)})
#' will assign PC1 to the X axis and PC2 to the Y axis. This parameter will be
#' ignored if X and Y are specified directly.
#' @param X directly supply the X coordinates, e.g. \code{pca_res$PC1}; if this
#' is used, \code{which_reddim} and \code{which_pcs} are ignored.
#' \code{X} should be of the same length (and order!) as \code{colData(object)}
#' (or \code{pData2.0(object)}).
#' @param Y directly supply the Y coordinates, e.g. \code{pca_res$PC2}.
#' Should be of the same length (and order!) as \code{colData(object)} (or
#' \code{pData2.0(object)}).
#' \code{which_reddim} and \code{which_pcs} will be ignored if \code{X} and 
#' \code{Y} are specified.
#' @param dim_red_type what type of dimensionality reduction was used to generate
#' the data (e.g. "PCA"). Will mostly be used for labelling purposes. Default: NULL
#' @param color_by entry of either \code{colData(object)}, \code{rowData(object)},
#' \code{rownames(object)}, which will be used to assign either discrete or
#' continous color schemes. Alternatively, you can supply a \code{list} with a
#' label stored in \code{list$title} and values stored in \code{list$result}
#' (\code{result} should be as long as \code{dim(object)[2]}. Default: \code{NULL}
#' @param shape_by default: \code{NULL}
#' @param size_by default: \code{NULL}
#' @param circle_by Used for encircling points of interest; default: \code{NULL}
#' @param exprs_values types of expression values that will be used for plotting if a gene name is being
#' specified in either \code{color_by}, \code{shape_by}, \code{size_by}
#' @param add_cell_info if you want to add columns for certain entries of
#' \code{colData(object)}, specify their names here. Default: NULL
#' 
#' @return DimRedPlot object: a simple S3-based object, i.e. a list with:
#' \enumerate{
#' \item \code{plot_data}: data.frame fit for ggplot2-based plotting
#' \item \code{label_data}: additional stuff that may be used for labelling the
#' x and y axes of the final plot
#' \item \code{factors}: names of the plotting variables, e.g. which values
#' were used for defining the coloring etc.
#' }
#' 
#' @seealso \code{\link{SingleCellExperiment}}, \code{\link{reducedDims}},
#' \code{\link{plt.DimRedPlot}}
#'
#' @import scater
#' 
#' @export
get_reducedDimPlot.sce = function(object,
    which_reddim = NULL,
    which_pcs = c(1,2),
    X = NULL, 
    Y = NULL,
    dim_red_type = NULL,
    color_by=NULL, 
    shape_by=NULL, 
    size_by=NULL,
    circle_by = NULL,
    exprs_values="counts",
    add_cell_info = NULL){
  
    if(is.null(which_reddim)  & (is.null(X) | is.null(Y)) ){
        stop("You must either provide the X and Y coordinates directly or specify the
            accessor for the reduced dimension coordinates via which_reddim.")
    }
  
    ## Define data.frame for plotting --------------------------------------------
  
    ## Either by using the cell names and the X and Y coordinates directly specified
    ## or by using the results stored within reducedDimension(object)
    if(!is.null(which_reddim)){
        use_reddim <- TRUE
    }
  
    if(!is.null(X) | !is.null(Y) ){
        if(any(is.null(X), is.null(Y))){
            # if only one of the axes is specified, we may resort to using which_reddim
            if(!is.null(which_reddim)){
                warning("You only supplied either X or Y. We'll be using which_reddim instead.")
            }else{
                stop("You must supply either which_reddim or both X and Y. Currently you're
                    only providing either X or Y.")
            }
        }else{
            use_reddim <- FALSE
        }
    }
  
    ## let's define the data.frame using the specific X and Y coordinates---------
    if(!use_reddim){
          xlen <- fx.return_axis_length(X)
          ylen <- fx.return_axis_length(Y)
          if(xlen != length(colnames(object)) | ylen != length(colnames(object))){
              stop("X and Y must be of the same length (and order!) as colnames(object).")
          }else{  
              ## make data.frame
              df_to_plot <- data.frame(x_axs = X, y_axs = Y, row.names = colnames(object))
      
              ## extract axes labels
              if( all(names(df_to_plot) == c("x_axs", "y_axs")) ){
                  axs_labs <- NULL
              }else{
                  axs_labs <- names(df_to_plot)
                  names(df_to_plot) <- c("x_axs", "y_axs")
              }
              rd <- NULL
          }
    }else{ ## Extracting coordinates from the reducedDims() slot -----------------
    
        if(!any(class(object) %in% "SingleCellExperiment")){
            stop("The object is not of class SingleCellExperiment, therefore we cannot retrieve 
                the reduced dimensions via reducedDims. Please provide X and Y directly.")
        }else{
            red_dim_all <- reducedDims(object)
        }
    
        ## check that there was actual data in that slot.................
        if(length(red_dim_all) == 0){
            stop("The object doesn't seem to have the reducedDim() slot filled.")
        }else{
            rd <- red_dim_all[[which_reddim]]
        }
    
        ## check the accessor....................
        if(is.null(rd)){
            stop("The accessor you provided for retrieving the coordinates from the
                reducedDims() list did not return any entries. Check that the name
                or index exists for reducedDims(object). Otherwise use X and Y to
                specify the vectors for the coordinates manually.")
        }
    
        ## finally, make the data.frame with coordinates for the components of interest...........
    
        if(is.null(which_pcs)){
            message("You didn't specify which components of the reduced dimension matrix you wanted to use.
                We're going to use the first two columns. If you want to change that, use which_pcs to
                specify the columns of interest, either via their names or their indexes.")
            which_pcs <- c(1,2)
        }else{
            ## make sure the specified coordinates are part of the reduced dimension object
            if( !is.numeric(which_pcs)){
                missing_pc <- which_pcs[ ! (which_pcs %in% colnames(rd)) ]
                if( length(missing_pc) > 0){
                    stop(paste0("The component(s) that you specified (",  paste(missing_pc, collapse=", "),") aren't part of the reduced dimension
                        object. That object contains the following possible entries: ", paste(colnames(rd), collapse = ", ")))
                }
            }
      
            if(length(which_pcs) > 2){
                warning(paste0("You indicated that you wanted to retrieve more than two components (", which_pcs, "). 
                   Since this function is geared towards x-y-plots, only the first two components will be used for the x and
                   y axis values, respectively."))
                which_pcs <- which_pcs[c(1:2)]
            }
        }
    
        ## subset the redDim data.frame
        df_to_plot <- data.frame(rd[, which_pcs, drop=FALSE])
    
        ## retrieve the original column names
        axs_labs <- names(df_to_plot)
    
        ## assign the standardized column names
        names(df_to_plot) <- c("x_axs","y_axs")
    }
  
    ## add meta-data to the X and Y coordinates, which will be used to assign color, shape etc. ------------
    if( !is.null(color_by) | !is.null(shape_by)  | !is.null(size_by) | !is.null(circle_by) ){
        df_to_plot <- fx.add_factors_to_redDim_df(object, 
            df_to_plot, 
            color_by, 
            shape_by,
            size_by, 
            circle_by,
            exprs_values)
    }
  
    ## return the data.frame with more info than what would normally be added based on color_by etc.
    ## this is useful if the data.frame is to be used in various ways without having to re-generate
    ## it all the time (e.g., extracting the values for individual genes can be time-consuming)
    if( !is.null(add_cell_info) ){
        df_to_plot <- fx.add_cellInfo_to_redDim_df(object, df_to_plot, add_cell_info)
    }
  
    ## gather data for labels in the plot etc.------------------------------------
    ld <- list()
    ld$x_lab <- axs_labs[1]
    ld$y_lab <- axs_labs[2]
    ld$variance_pct <- fx.get_percent_var(rd,  pcs_select = which_pcs )
    ld$exprs_val_type <- exprs_values
    ld$dim_red_type <- dim_red_type
  
    ## remember which factors were used for what type of plotting variable--------
    fc <- list()
    fc$color_by <- color_by
    fc$size_by <- size_by
    fc$shape_by <- shape_by
    fc$circle_by <- circle_by
  
    ## return an S3 object--------------------------------------------------------
    DRP <- list(plot_data = df_to_plot,
        label_data = ld,
        factors = fc)
    class(DRP) <- "DimRedPlot"
  
    return(DRP)
}
