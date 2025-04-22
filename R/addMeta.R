#' Add a metadata to be included in the shiny app
#'
#' Add a metadata to be included in the shiny app.
#'
#' @param scConf shinycell config data.table
#' @param meta.to.add metadata to add from the single-cell metadata.  
#'   Must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: column names in \code{seu@meta.data} 
#'       i.e. \code{colnames(seu@meta.data)}
#'     \item{h5ad files}: column names in \code{h5ad.obs} 
#'       i.e. \code{h5ad.obs.columns.values} 
#'   }
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param maxLevels maximum number of levels allowed for categorical metadata.
#'   Metadata with nlevels > maxLevels will throw up an error message
#' 
#' @return updated shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table reticulate hdf5r
#'
#' @examples
#' scConf = addMeta(scConf, c("orig.ident"), seu)
#'
#' @export
addMeta <- function(scConf, meta.to.add, obj, maxLevels = 50){
  # Extract corresponding metadata
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    objMeta = obj@meta.data
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    ad <- import("anndata", convert = FALSE)
    inpH5 = ad$read_h5ad(obj)
    objMeta = data.frame(py_to_r(inpH5$obs$values))
    rownames(objMeta) = py_to_r(inpH5$obs_names$values)
    colnames(objMeta) = py_to_r(inpH5$obs$columns$values)
    for(i in colnames(objMeta)){
      objMeta[[i]] = unlist(objMeta[[i]])   # unlist and refactor
      if(as.character(inpH5$obs[i]$dtype) == "category"){
        objMeta[[i]] = factor(objMeta[[i]], levels = 
                                py_to_r(inpH5$obs[i]$cat$categories$values))
      }
    } 
    
  } else {
    stop("Only Seurat objects or h5ad file paths are accepted!")
  }
  
  # Check if meta.to.add exist in obj metadata
  if(!all(meta.to.add %in% colnames(objMeta))){
    stop("meta.to.add not found in single-cell data object!")
  }
  
  # Start adding meta.data
  colCateg = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
               "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  for(iMeta in meta.to.add){
    tmpConf = data.table(ID = iMeta, UI = iMeta, fID = NA, fUI = NA, 
                         fCL = NA, fRow = NA, default = 0, grp = FALSE)
    
    # Convert to factors if metadata contains characters
    if(is.character(objMeta[[iMeta]])){
      objMeta[[iMeta]] = factor(objMeta[[iMeta]])
    }
    
    # Additional preprocessing for categorical metadata
    nLevels = nlevels(objMeta[[iMeta]])
    if(nLevels > maxLevels){
      stop(paste0(iMeta, " has exceeded the maximum number of levels allowed!"))
    }
    if(nLevels >= 2){
      tmpConf$fID = paste0(levels(objMeta[[iMeta]]), collapse="|")
      tmpConf$fUI = tmpConf$fID
      tmpConf$fCL = paste0(colorRampPalette(colCateg)(nLevels), collapse="|")
      tmpConf$fRow = ceiling(nLevels / 4)
      tmpConf$grp = TRUE
    } else if(nLevels == 1){
      tmpConf$fID = levels(objMeta[[iMeta]])
      tmpConf$fUI = tmpConf$fID
      tmpConf$fCL = "black"
      tmpConf$fRow = 1
    }
    scConf = rbindlist(list(scConf, tmpConf))
  }
  return(scConf)
}


