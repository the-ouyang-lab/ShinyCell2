#' Make h5 object from Seurat assay data
#'
#' Make h5 object from Seurat assay data
#'
#' @param obj input Seurat (v3+) object
#' @param sc1meta data.table of cell metadata
#' @param filename filename of output h5 file
#' @param gex.assay assay in Seurat object to use
#' @param gex.slot slot in single-cell assay to use
#' @param chunkSize number of genes written to h5file at any one time
#'
#' @return h5 object
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate
#'
#' @export
makeH5fromSeurat <- function(obj, sc1meta, filename, 
                             gex.assay, gex.slot, chunkSize){
  # Create h5 file and get ready
  if(class(obj@assays[[gex.assay]]) == "Assay5"){
    gex.matdim = dim(obj@assays[[gex.assay]]@layers[[gex.slot]])
  }else{
    gex.matdim = dim(slot(obj@assays[[gex.assay]], gex.slot))
  }
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",  dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1,gex.matdim[2]))
  chk = chunkSize
  while(chk > (gex.matdim[1]-8)){
    chk = floor(chk / 2)     # Account for cases where nGene < chunkSize
  } 
  
  # Start writing to file
  nChunk = floor((gex.matdim[1]-8)/chk)
  if(class(obj@assays[[gex.assay]]) == "Assay5"){
    # Seurat v5
    for(i in 1:nChunk){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        obj[[gex.assay]][gex.slot][
          ((i-1)*chk+1):(i*chk), sc1meta$cellID])
    }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        obj[[gex.assay]][gex.slot][
          (i*chk+1):gex.matdim[1], sc1meta$cellID])
      gex.rownm = rownames(obj[[gex.assay]][gex.slot])
      
  } else {
    for(i in 1:nChunk){
      sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
        slot(obj@assays[[gex.assay]], gex.slot)[
          ((i-1)*chk+1):(i*chk), sc1meta$cellID])
    }
      sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
        slot(obj@assays[[gex.assay]], gex.slot)[
          (i*chk+1):gex.matdim[1], sc1meta$cellID])
      gex.rownm = rownames(slot(obj@assays[[gex.assay]], gex.slot))
  }
  sc1gexpr$close_all()
  # Output gene names
  return(gex.rownm)
}


