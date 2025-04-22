#' Create a shinycell config data.table
#'
#' Create a shinycell config data.table containing (i) the single-cell 
#' metadata to display on the Shiny app, (ii) ordering of factors / 
#' categories of categorical metadata and (iii) colour palettes associated 
#' with each metadata.
#'
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param meta.to.include columns to include from the single-cell metadata. 
#'   Default is \code{NA}, which is to use all columns. Users can specify 
#'   the columns to include, which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: column names in \code{seu@meta.data} 
#'       i.e. \code{colnames(seu@meta.data)}
#'     \item{h5ad files}: column names in \code{h5ad.obs} 
#'       i.e. \code{h5ad.obs.columns.values} 
#'   }
#' @param legendCols maximum number of columns allowed when displaying the 
#'   legends of categorical metadata
#' @param maxLevels maximum number of levels allowed for categorical metadata.
#'   Metadata with nlevels > maxLevels will be discarded automatically
#'
#' @return shinycell config data.table
#'
#' @author John F. Ouyang
#'
#' @import data.table reticulate hdf5r
#'
#' @examples
#' scConf = createConfig(obj)
#'
#' @export
createConfig <- function(obj, meta.to.include = NA, legendCols = 4,
                         maxLevels = 50){
  # Extract corresponding metadata
  drExist = TRUE
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    objMeta = obj@meta.data
    if(length(names(obj@reductions)) == 0){drExist = FALSE}
    
  } else if(class(obj)[1] == "ArchRProject"){
    # ArchR Object
    objMeta = as.data.frame(obj@cellColData)
    if(length(names(obj@embeddings@listData)) == 0){drExist = FALSE}
    
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
    if(length(py_to_r(inpH5$obsm_keys())) == 0){drExist = FALSE}

  } else {
    stop("Only Seurat/ArchR objects or h5ad file paths are accepted!")
  }
  if(!drExist){
    stop(paste0("ShinyCell did not detect any dimension reduction data \n", 
                "       e.g. umap / tsne. Has any analysis been performed?"))
  }
  
  # Checks and get list of metadata to include
  if(is.na(meta.to.include[1])){meta.to.include = colnames(objMeta)}
  if(length(meta.to.include) < 2){stop("At least 2 metadata is required!")}
  if(!all(meta.to.include %in% colnames(objMeta))){stop("Some metadata are missing!")}
  
  # Start making config data.table
  colCateg = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C",
               "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  scConf = data.table()
  for(iMeta in meta.to.include){
    tmpConf = data.table(ID = iMeta, UI = iMeta, fID = NA, fUI = NA, 
                         fCL = NA, fRow = NA, default = 0, grp = FALSE)
    
    # Convert to factors if metadata contains characters
    if(is.character(objMeta[[iMeta]])){
      objMeta[[iMeta]] = factor(objMeta[[iMeta]])
    }
    
    # Additional preprocessing for categorical metadata
    nLevels = nlevels(objMeta[[iMeta]])
    if(nLevels <= maxLevels){
      if(nLevels >= 2){
        tmpConf$fID = paste0(levels(objMeta[[iMeta]]), collapse="|")
        tmpConf$fUI = tmpConf$fID
        tmpConf$fCL = paste0(colorRampPalette(colCateg)(nLevels), collapse="|")
        tmpConf$fRow = ceiling(nLevels / legendCols)
        tmpConf$grp = TRUE
      } else if(nLevels == 1){
        tmpConf$fID = levels(objMeta[[iMeta]])
        tmpConf$fUI = tmpConf$fID
        tmpConf$fCL = "black"
        tmpConf$fRow = 1
      }
      scConf = rbindlist(list(scConf, tmpConf))
    } else {
      warning(paste0(iMeta, " is excluded as it has >", maxLevels, " levels!"))
    }
  }
  
  # Set defaults
  def1 = grep("ident|library", scConf$ID, ignore.case = TRUE)[1]
  def2 = grep("clust", scConf$ID, ignore.case = TRUE)
  def2 = setdiff(def2, def1)[1]
  if(is.na(def1)){def1 = setdiff(c(1,2), def2)[1]}
  if(is.na(def2)){def2 = setdiff(c(1,2), def1)[1]}
  scConf[def1]$default = 1
  scConf[def2]$default = 2
  
  # STOP if there is no single multi-level covariate
  if(nrow(scConf[grp == TRUE]) == 0){
    stop(paste0("ShinyCell did not detect any multi-group cell metadata \n", 
                "       e.g. library / cluster. Has any analysis been performed?"))
  }
  
  return(scConf)
}


