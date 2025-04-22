#' Generate data files required for shiny app (scRNA type data)
#'
#' Generate data files required for shiny app, specifically scRNA-seq data
#' Six files will be generated, namely 
#' (i) the shinycell config \code{prefix_conf.rds}, 
#' (ii) the single-cell metadata \code{prefix_meta.rds}, 
#' (iii) the single-cell assays \code{prefix_assay_X.h5}, 
#' (iv) the feature mapping object config \code{prefix_gene.rds}, 
#' (v) the dimension reduction embeddings \code{prefix_dimr.rds} and 
#' (vi) the defaults for the Shiny app \code{prefix_def.rds}. 
#' A prefix is specified for each set of files to allow for multiple 
#' single-cell datasets in a single Shiny app.
#'
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param scConf shinycell config data.table
#' @param gex.assay assay(s) in single-cell data object to use. Multiple assays
#'   can now be incorporated and all assays are used by default (with the first 
#'   assay being the default assay), which must match one of the following:
#'   \itemize{
#'     \item{Seurat objects}: "RNA" or "integrated" assay, 
#'       default is "RNA"
#'     \item{h5ad files}: "X" or any assay in "layers",
#'       default is "X"
#'   }
#' @param gex.slot slot in single-cell assay to plot. This is only used 
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
#' makeShinyFilesGEX(seu, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/",
#'                   default.gene1 = "POU5F1", default.gene2 = "APOA1",
#'                   default.multigene = c("POU5F1","APOA1","GPRC5A","TBXT","ISL1"),
#'                   default.dimred = "umap")
#'
#' @export
makeShinyFilesGEX <- function(
  obj, scConf, gex.assay = NA, gex.slot = "data", dimred.to.use = NA,
  shiny.prefix = "sc1", shiny.dir = "shinyApp/",
  default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
  default.dimred = NA, chunkSize = 500){
  ### Preprocessing and checks
  # Generate defaults for gex.assay / gex.slot
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    if(is.na(gex.assay[1])){
      gex.assay = names(obj@assays)
      if (requireNamespace("SeuratObject", quietly = TRUE)){
        gex.assay = c(SeuratObject::DefaultAssay(obj), 
                      setdiff(gex.assay, SeuratObject::DefaultAssay(obj)))
      }
      gex.assay = setdiff(gex.assay, "peaks")
    }
    if(class(obj@assays[[gex.assay[1]]]) == "Assay5"){
      # Seurat v5: check if layers are joined
      if(!(gex.slot[1] %in% names(obj@assays[[gex.assay[1]]]@layers))){
        stop(paste0("gex.slot not found in gex.assay. ",
                    "Are layers joined? run obj <- JoinLayers(obj)"))
      }
    }
    defGenes = Seurat::VariableFeatures(obj)[1:10]
    if(is.na(defGenes[1])){
      warning(paste0("Variable genes for seurat object not found! Have you ",
                     "ran `FindVariableFeatures` or `SCTransform`?"))
      defGenes = rownames(obj)[1:10]
    }
    if(is.na(dimred.to.use[1])){
      dimred.to.use = setdiff(names(obj@reductions), c("pca","ref.pca","lsi"))
      dimred.to.use = c(DefaultDimReduc(obj), setdiff(dimred.to.use, DefaultDimReduc(obj)))
    }
    sc1meta = data.table(cellID = rownames(obj@meta.data), obj@meta.data)
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    if(is.na(gex.assay[1])){gex.assay = "X"}
    # Can just check X since inpH5$layers should have same dimensions
    ad <- import("anndata", convert = FALSE)
    sp <- import('scipy.sparse', convert = FALSE)
    inpH5 = ad$read_h5ad(obj)
    defGenes = py_to_r(inpH5$var_names$values)[1:10]
    if(is.na(dimred.to.use[1])){
      dimred.to.use = setdiff(py_to_r(inpH5$obsm_keys()), c("X_pca"))}
    sc1meta = data.table(cellID = py_to_r(inpH5$obs_names$values))
    sc1meta = cbind(sc1meta, data.table(py_to_r(inpH5$obs$values)))
    colnames(sc1meta) = c("cellID", py_to_r(inpH5$obs$columns$values))
    for(i in colnames(sc1meta)[-1]){
      sc1meta[[i]] = unlist(sc1meta[[i]])   # unlist and refactor
      if(as.character(inpH5$obs[i]$dtype) == "category"){
        sc1meta[[i]] = factor(sc1meta[[i]], levels = 
                                py_to_r(inpH5$obs[i]$cat$categories$values))
      }
    } 

  } else {
    stop("Only Seurat objects or h5ad file paths are accepted!")
  }
  # Check default.gene1 / default.gene2 / default.multigene
  if(is.na(default.gene1[[1]])){default.gene1 = defGenes[1]}
  if(is.na(default.gene2[[1]])){default.gene2 = defGenes[2]}
  if(is.na(default.multigene[[1]])){default.multigene = defGenes}
  default.dimred = default.dimred[1]
  if(is.na(default.dimred)){default.dimred = dimred.to.use[1]}
  

  
  ### Actual object generation
  # Make prefix_conf.rds / prefix_meta.rds
  sc1conf = scConf
  sc1meta = sc1meta[, c("cellID", as.character(sc1conf$ID)), with = FALSE]
  # Factor metadata again
  for(i in as.character(sc1conf[!is.na(fID)]$ID)){
    sc1meta[[i]] = factor(sc1meta[[i]],
                          levels = strsplit(sc1conf[ID == i]$fID, "\\|")[[1]])
    levels(sc1meta[[i]]) = strsplit(sc1conf[ID == i]$fUI, "\\|")[[1]]
    sc1conf[ID == i]$fID = sc1conf[ID == i]$fUI
  }
  sc1conf$ID = as.character(sc1conf$ID)     # Remove levels
  
  # Make prefix_dimr.rds
  sc1dimr = list()
  if(class(obj)[1] == "Seurat"){
    # Seurat Object
    for(iDR in dimred.to.use){
      drMat = obj@reductions[[iDR]]@cell.embeddings
      if(ncol(drMat) > 2){drMat = drMat[, 1:2]}  # Take first two comps only
      drMat = drMat[sc1meta$cellID, ]            # Ensure ordering
      sc1dimr[[iDR]] = drMat
    }
    
  } else if (tolower(tools::file_ext(obj)) == "h5ad"){
    # h5ad file
    for(iDR in dimred.to.use){
      drMat = py_to_r(inpH5$obsm[iDR])
      tmpName = gsub("pca", "pc", gsub("X_", "", iDR))
      tmpName = paste0(tmpName, "_", 1:ncol(drMat))
      colnames(drMat) = tmpName
      if(ncol(drMat) > 2){drMat = drMat[, 1:2]}  # Take first two comps only
      drMat = drMat[sc1meta$cellID, ]            # Ensure ordering
      sc1dimr[[iDR]] = drMat
    }
    
  }
  
  # Make prefix_assay_X.h5 / prefix_gene.rds
  sc1gene = list()
  if(!dir.exists(shiny.dir)){dir.create(shiny.dir)}
  for(iAssay in gex.assay){
    filename = paste0(shiny.dir, "/", shiny.prefix, "assay_", iAssay, ".h5")
    if(class(obj)[1] == "Seurat"){
      tmpG = makeH5fromSeurat(obj, sc1meta, filename = filename,
                              gex.assay = iAssay, gex.slot = gex.slot[1], 
                              chunkSize = chunkSize)
    } else if (tolower(tools::file_ext(obj)) == "h5ad"){
      tmpG = makeH5fromAnndata(obj, sc1meta, filename = filename, 
                               gex.assay = iAssay, chunkSize = chunkSize)
    }
    tmpOut = seq_along(tmpG); names(tmpOut) = tmpG
    tmpOut = tmpOut[order(names(tmpOut))]
    tmpOut = tmpOut[order(nchar(names(tmpOut)))]
    sc1gene[[iAssay]] = tmpOut
  }
  # Check (again) default.gene1 / default.gene2 / default.multigene
  defGene1 = list()
  defGene2 = list()
  defGeneM = list()
  for(iAssay in gex.assay){
    tmp0 = names(sc1gene[[iAssay]][which(sc1gene[[iAssay]] %in% 1:10)])
    # default.gene1
    if(is.list(default.gene1) & (length(default.gene1) >= length(gex.assay))){
      tmp = default.gene1[[iAssay == gex.assay]][1]
    } else {
      tmp = default.gene1[[1]][1]
    }
    if(tmp %in% names(sc1gene[[iAssay]])){
      defGene1[[iAssay]] = tmp
    } else {
      defGene1[[iAssay]] = tmp0[1]
    }
    # default.gene2
    if(is.list(default.gene2) & (length(default.gene2) >= length(gex.assay))){
      tmp = default.gene2[[iAssay == gex.assay]][1]
    } else {
      tmp = default.gene2[[1]][1]
    }
    if(tmp %in% names(sc1gene[[iAssay]])){
      defGene2[[iAssay]] = tmp
    } else {
      defGene2[[iAssay]] = tmp0[2]
    }
    # default.multigene
    if(is.list(default.multigene) & (length(default.multigene) >= length(gex.assay))){
      tmp = default.multigene[[iAssay == gex.assay]]
    } else {
      tmp = unique(unlist(default.multigene))
    }
    if(all(tmp %in% names(sc1gene[[iAssay]]))){
      defGeneM[[iAssay]] = tmp
    } else {
      defGeneM[[iAssay]] = tmp0
    }
  }
  
  # Make prefix_def.rds
  # Note that we stored the display name here
  sc1def = list()
  sc1def$meta1 = sc1conf[default == 1]$UI   # Use display name
  sc1def$meta2 = sc1conf[default == 2]$UI   # Use display name 
  sc1def$gene1 = defGene1         # Actual == Display name
  sc1def$gene2 = defGene2         # Actual == Display name
  sc1def$genes = defGeneM         # Actual == Display name
  sc1def$dimrd = dimred.to.use    # First one is default    
  sc1def$assay = gex.assay        # First one is default
  tmp = nrow(sc1conf[default != 0 & grp == TRUE])
  if(tmp == 2){
    sc1def$grp1 = sc1def$meta1
    sc1def$grp2 = sc1def$meta2
  } else if(tmp == 1){
    sc1def$grp1 = sc1conf[default != 0 & grp == TRUE]$UI
    if(nrow(sc1conf[default == 0 & grp == TRUE]) == 0){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    }
  } else {
    sc1def$grp1 = sc1conf[default == 0 & grp == TRUE]$UI[1]
    if(nrow(sc1conf[default == 0 & grp == TRUE]) < 2){
      sc1def$grp2 = sc1def$grp1
    } else {
      sc1def$grp2 = sc1conf[default == 0 & grp == TRUE]$UI[2]
    }
  }
  sc1conf = sc1conf[, -c("fUI", "default"), with = FALSE]
  
  
  
  ### Saving objects
  saveRDS(sc1conf, file = paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))
  saveRDS(sc1meta, file = paste0(shiny.dir, "/", shiny.prefix, "meta.rds"))
  saveRDS(sc1gene, file = paste0(shiny.dir, "/", shiny.prefix, "gene.rds"))
  saveRDS(sc1dimr, file = paste0(shiny.dir, "/", shiny.prefix, "dimr.rds"))
  saveRDS(sc1def,  file = paste0(shiny.dir, "/", shiny.prefix, "def.rds"))
  return(paste0(shiny.prefix, " files generated!"))
}


