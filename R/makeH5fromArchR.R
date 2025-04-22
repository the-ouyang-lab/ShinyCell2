#' Make h5 object from ArchR data
#'
#' Make h5 object from ArchR data
#'
#' @param obj input ArchR object
#' @param sc1meta data.table of cell metadata
#' @param filename filename of output h5 file
#' @param assay assay in ArchR object to use
#'
#' @return h5 object
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate 
#'
#' @export
makeH5fromArchR <- function(obj, sc1meta, filename, assay){
  # Extract peakset and generate chrList 
  chrList <- levels(seqnames(obj@peakSet))
  # Create h5 file and get ready
  if(assay == "PeakMatrix"){
    gex.rownm <- obj@peakSet
    gex.rownm <- paste0(seqnames(gex.rownm),":", 
                        start(gex.rownm),"-",end(gex.rownm))
  } else {
    gex.rownm <- getFeatures(ArchRProj = obj, useMatrix = assay)
  }
  if(assay == "MotifMatrix"){gex.rownm <- gex.rownm[grep("^z:", gex.rownm)]}
  gex.matdim = c(length(gex.rownm), nrow(sc1meta))
  sc1gexpr <- H5File$new(filename, mode = "w")
  sc1gexpr.grp <- sc1gexpr$create_group("grp")
  sc1gexpr.grp.data <- sc1gexpr.grp$create_dataset(
    "data",  dtype = h5types$H5T_NATIVE_FLOAT,
    space = H5S$new("simple", dims = gex.matdim, maxdims = gex.matdim),
    chunk_dims = c(1,gex.matdim[2]))

  # Start writing to file
  if(assay == "MotifMatrix"){
    tmp <- getMatrixFromProject(ArchRProj = obj, useMatrix = "MotifMatrix")
    tmp <- as.matrix(tmp@assays@data@listData[["z"]])
    sc1gexpr.grp.data[,] <- tmp[, sc1meta$cellID]
    gex.rownm <- rownames(tmp)
  } else {
    gex.rownm = c()
    # Extract by chr for peak/gene matrices
    for (iChr in chrList) {
      # Generate summarizedexperiment and collect rownames
      tmp <- getMatrixFromProject(ArchRProj = obj,
                                  useMatrix = assay,
                                  useSeqnames = iChr)
      if(assay == "PeakMatrix"){
        tmp.rownm <- paste0(seqnames(tmp@rowRanges),":", 
                            start(tmp@rowRanges),"-",end(tmp@rowRanges))
      } else {
        tmp.rownm <- rowData(tmp)$name
      }
      gex.rownm <- c(gex.rownm, tmp.rownm)
      # Subset out peakmatrix and add to h5file
      tmp <- as.matrix(tmp@assays@data@listData[[assay]])
      sc1gexpr.grp.data[(length(gex.rownm)-length(tmp.rownm)+1):(
        length(gex.rownm)),] <- tmp[, sc1meta$cellID]
    }
  }
  sc1gexpr$close_all()
  
  # Output gene names
  return(gex.rownm)
}


