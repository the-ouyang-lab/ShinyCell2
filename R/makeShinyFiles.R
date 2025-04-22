#' Generate data files required for shiny app
#'
#' Wrapper function to generate data files required for shiny app. 
#' A prefix is specified for each set of files to allow for multiple 
#' single-cell datasets in a single Shiny app.
#'
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param scConf shinycell config data.table
#' @param bigWigGroup categorical group in scATAC datasets to group cells by for 
#'   the generation of bigWig files for track plot. Default is NA which does 
#'   not generate any bigWig files.
#' @param assay assay(s) in single-cell data object to use. Multiple assays
#'   can now be incorporated and all assays are used by default (with the first 
#'   assay being the default assay), which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay, 
#'       default is "RNA"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'   }
#' @param assay.slot slot in single-cell assay to plot. This is only used 
#'   for Seurat objects (v3+). Default is to use the "data" slot
#' @param dimred.to.use specify the dimension reduction to use. Default is to 
#'   use all except PCA 
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default gene to show, which be present 
#'   in the default assay for Seurat or X layer in scanpy h5ad
#' @param default.gene2 specify secondary default gene to show, which be present 
#'   in the default assay for Seurat or X layer in scanpy h5ad
#' @param default.multigene character vector specifying default genes to 
#'   show in bubbleplot / heatmap, which be present 
#'   in the default assay for Seurat or X layer in scanpy h5ad
#' @param default.dimred character vector specifying the two default dimension 
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param chunkSize number of genes written to h5file at any one time. Lower 
#'   this number to reduce memory consumption. Should not be less than 10
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate hdf5r
#'
#' @examples
#' makeShinyFiles(seu, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/",
#'                default.gene1 = "POU5F1", default.gene2 = "APOA1",
#'                default.multigene = c("POU5F1","APOA1","GPRC5A","TBXT","ISL1"),
#'                default.dimred = "umap")
#'
#' @export
makeShinyFiles <- function(
    obj, scConf, bigWigGroup = NA, 
    assay = NA, assay.slot = "data", dimred.to.use = NA,
    shiny.prefix = "sc1", shiny.dir = "shinyApp/",
    default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
    default.dimred = NA, chunkSize = 500, ...){
  ### Check object class
  if(class(obj)[1] == "Seurat"){
    # Seurat object
    makeShinyFilesGEX(
      obj, scConf, assay, assay.slot, dimred.to.use,
      shiny.prefix, shiny.dir,
      default.gene1, default.gene2, default.multigene, 
      default.dimred, chunkSize)
    # Seurat object with spatial data
    if(.hasSlot(obj, "images")){
      if(length(obj@images) > 0){
        makeShinyFilesSpatial(
          obj, scConf, shiny.prefix, shiny.dir)
      }
    }
    # Signac object with fragment data
    if("peaks" %in% names(obj@assays)){
      if(!is.na(bigWigGroup[1])){
        makeShinyFilesATACsignac(
          obj, scConf, bigWigGroup, shiny.prefix, shiny.dir)
      }
    }

  } else if (class(obj)[1] == "ArchRProject"){
    # ArchR object
    makeShinyFilesATACarchr(
      obj, scConf, bigWigGroup, assay, dimred.to.use,
      shiny.prefix, shiny.dir,
      default.gene1, default.gene2, default.multigene, 
      default.dimred, chunkSize, ...)
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    makeShinyFilesGEX(
      obj, scConf, assay, assay.slot, dimred.to.use,
      shiny.prefix, shiny.dir,
      default.gene1, default.gene2, default.multigene, 
      default.dimred, chunkSize)
    
  }
}


