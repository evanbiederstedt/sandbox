
#' @name my_gene
#' @title Map between ENSEMBL IDs and gene symbols
#' @description
#' 
#' Convenience function for quickly extracting the gene symbol for
#' an ENSEMBL ID of interest or vice versa. This is needed often for the Lo data
#' set because the row names of the SCE object are ENSEMBL IDs.
#' 
#' @param goi vector of the IDs of the gene(s) of interest, e.g. "Ins2"
#' @param gn_map object that contains the mapping between gene_symbol and ensembl_id,
#' either an SCE object (from which \code{rowData()} will be used) or a data.table
#' with "gene_symbol" and "ensembl_id"
#' 
#' @return either the gene_symbols (if ENSEMBL IDs were given) or the ENSEMBL IDs
#' (if gene_symbols were given) of the genes of interest
#' 
#' @examples \dontrun{
#' gns <- rowData(sce.filt)[c(1:2)] %>% as.data.table
#' my_gene("Ins2", gns)
#' # [1] "ENSMUSG00000000215"
#' }
#' 
#' @import data.table
#' 
#' @export
#' 
my_gene = function(goi, gn_map){
  
    if(any(class(gn_map) == "SingleCellExperiment")){
        rd <- rowData(gn_map)
        ABCutilities::check_columns(c("gene_symbol","ensembl_id"), rd, "gn_map", "my_gene")
        gnmp <- as.data.table(rd[c("gene_symbol","ensembl_id")])
    }else{
        if(any(class(gn_map) == "data.table")){
            ABCutilities::check_columns(c("gene_symbol","ensembl_id"), gn_map, "gn_map", "my_gene")
            gnmp <- gn_map
        }else{
            stop("The gene map is neither a SCE object nor a data.table.")
        }
    }
  
    if(all(grepl("^ENSMUS", goi))){
        out <- unique(gnmp[ensembl_id %in% goi]$gene_symbol)
    }else if(all(!grepl("^ENSMUS", goi))){
        out <- unique(gnmp[gene_symbol %in% goi]$ensembl_id)
    }else{
        stop("All genes of interest should be of the same ID type, either ENSEMBL IDs or gene_symbols only.")
    }

    return(out)
    
}
