#' Make h5 object from Anndata assay
#'
#' Make h5 object from Anndata assay
#'
#' @param obj input file path for h5ad file
#' @param sc1meta data.table of cell metadata
#' @param filename filename of output h5 file
#' @param gex.assay assay in anndata object to use
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
  ad <- import("anndata", convert = FALSE)
  sp <- import('scipy.sparse', convert = FALSE)
  inpH5 = ad$read_h5ad(obj)
  gex.matdim = rev(unlist(py_to_r(inpH5$X$shape)))  
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
  if(gex.assay == "X"){
    scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$X)))
  } else {
    scGEX = Matrix::t(py_to_r(sp$csc_matrix(inpH5$layers[[gex.assay]])))
  }
  for(i in 1:nChunk){
    sc1gexpr.grp.data[((i-1)*chk+1):(i*chk), ] <- as.matrix(
      scGEX[((i-1)*chk+1):(i*chk), sc1meta$cellID])
  }
  sc1gexpr.grp.data[(i*chk+1):gex.matdim[1], ] <- as.matrix(
    scGEX[(i*chk+1):gex.matdim[1], sc1meta$cellID])
  sc1gexpr$close_all()
  
  # Output gene names
  gex.rownm = py_to_r(inpH5$var_names$values)
  return(gex.rownm)
}


