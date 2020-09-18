#' This function takes in a Seurat 3.0 object and returns a named list with 2
#' objects formated to be loaded in the ShinyApp:
#'
#'      1. Metadata + coordinates of the 2D embeding.
#'      2. Expression matrix of the selected slot.
#'
#' ShinyApp: https://singlecellgenomics-cnag-crg.shinyapps.io/Annotation
#'
#' @param object Object of class Seurat. Mandatory.
#' @param assay Character string. Assay within the Seurat object from which to extract the expression matrix. Default: active assay.
#' @param reduction Character string. Dimensionality reduction from which to extract the 2D coordinates. Default: umap.
#' @param slot Character string. Slot containing the expression matrix. Default: data.
#' @param asfactors Character vector. Metadata columns to be converted to factors. Default: NULL.
#' @param save Logical value. Save metadata and expression matrix objects as RDS files. Default: TRUE.
#' @param path Character string. Path to save output files if 'save = TRUE'. Default: working directory.
#'
#' @return Named list containing the joint metadata + 2D embeding and the expression matrix.
#'
#' @examples
#' seurat2shiny( object = seurat_object, asfactors = c("plate", "replicate") )
#'
#' shiny_list = seurat2shiny(object = seurat_object)
#'
#' @export


seurat2shiny = function(
    object                         ,
    assay     = object@active.assay,
    reduction = "umap"             ,
    slot      = "data"             ,
    asfactors = NULL               ,
    save      = TRUE               ,
    path      = "."                  # path = getwd()
) {
    suppressMessages( library(Seurat) );

    # Input check.
    if ( ! is(object, "Seurat") )
        stop("'object' is not a Seurat object.");

    if ( ! assay %in% Seurat::Assays(object) )
        stop("'assay' not in the Seurat object's available assays.");

    if ( ! reduction %in% names(object@reductions) )
        stop("'reduction' not in the Seurat object's available reductions.");

    if ( ! slot %in% c("counts", "data", "scale.data") )
        stop("'slot' not in the Seurat object's available slots.");

    # Extract 2D coordinates.
    embeds = as.data.frame(object@reductions[[reduction]]@cell.embeddings);
    names(embeds) = c("coord_x", "coord_y");

    # Join metadata with coordinates.
    metadata = object@meta.data;

    for (col in asfactors) {
        metadata[[col]] = as.factor(metadata[[col]]);
    };

    metadata = merge(x = metadata, y = embeds, by = "row.names");
    names(metadata)[1] = "barcode"; # names(metadata)[names(metadata) == "Row.names"] = "barcode";

    # Extract expression data.
    # expression = as.matrix( Seurat::GetAssayData(object = object, slot = slot, assay = assay) );
    expression = Seurat::GetAssayData(object = object, slot = slot, assay = assay);

    if ( ! identical( as.character(metadata$barcode), colnames(expression) ) )
        warning("Cells in metadata and expression matrix do not match.");

    if (save) {
        saveRDS( object = metadata  , file = paste0(path, "/metadata.rds"  ) );
        saveRDS( object = expression, file = paste0(path, "/expression.rds") );
    };

    invisible(
        list(metadata = metadata, expression = expression)
    );
};
