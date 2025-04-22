ExportGroupBW  <- function(
    object,
    assay = NULL,
    group.by = NULL,
    idents = NULL,
    normMethod = "RC",
    tileSize = 100,
    minCells = 5,
    cutoff = NULL,
    chromosome = NULL,
    outdir = NULL,
    verbose=TRUE
) {
  # Check if temporary directory exist
  if (!dir.exists(outdir)){
    dir.create(outdir)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  assay <- SetIfNull(x = assay, y = DefaultAssay(object = object))
  DefaultAssay(object = object) <- assay
  group.by <- SetIfNull(x = group.by, y = 'ident')
  Idents(object = object) <- group.by
  idents <- SetIfNull(x = idents, y = levels(x = object))
  GroupsNames <- names(x = table(object[[group.by]])[table(object[[group.by]]) > minCells])
  GroupsNames <- GroupsNames[GroupsNames %in% idents]
  # Check if output files already exist
  lapply(X = GroupsNames, FUN = function(x) {
    fn <- paste0(outdir, .Platform$file.sep, x, ".bed")
    if (file.exists(fn)) {
      message(sprintf("The group \"%s\" is already present in the destination folder and will be overwritten !",x))
      file.remove(fn)
    }
  })      
  # Splitting fragments file for each idents in group.by
  SplitFragments(
    object = object,
    assay = assay,
    group.by = group.by,
    idents = idents,
    outdir = outdir,
    file.suffix = "",
    append = TRUE,
    buffer_length = 256L,
    verbose = verbose
  )
  # Column to normalized by
  if(!is.null(x = normMethod)) {
    if (tolower(x = normMethod) %in% c('rc', 'ncells', 'none')){
      normBy <- normMethod
    } else{
      normBy <- object[[normMethod, drop = FALSE]]
    }
  }
  # Get chromosome information
  if(!is.null(x = chromosome)){
    seqlevels(object) <- chromosome
  }
  availableChr <- names(x = seqlengths(object))
  chromLengths <- seqlengths(object)
  chromSizes <- GRanges(
    seqnames = availableChr,
    ranges = IRanges(
      start = rep(1, length(x = availableChr)),
      end = as.numeric(x = chromLengths)
    )
  )
  
  if (verbose) {
    message("Creating tiles")
  }
  # Create tiles for each chromosome, from GenomicRanges
  tiles <- unlist(
    x = slidingWindows(x = chromSizes, width = tileSize, step = tileSize)
  )
  if (verbose) {
    message("Creating bigwig files at ", outdir)
  }
  # Run the creation of bigwig for each cellgroups
  if (nbrOfWorkers() > 1) { 
    mylapply <- future_lapply 
  } else { 
    mylapply <- lapply
  }
  covFiles <- mylapply(
    GroupsNames,
    FUN = CreateBWGroup,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
  )
  return(covFiles)
}


CreateBWGroup <- function(
    groupNamei,
    availableChr,
    chromLengths,
    tiles,
    normBy,
    tileSize,
    normMethod,
    cutoff,
    outdir
) {
  if (!requireNamespace("rtracklayer", quietly = TRUE)) { 
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/") 
    return(NULL) 
  }
  normMethod <- tolower(x = normMethod)
  # Read the fragments file associated to the group
  fragi <- rtracklayer::import(
    paste0(outdir, .Platform$file.sep, groupNamei, ".bed"), format = "bed"
  )
  cellGroupi <- unique(x = fragi$name)
  # Open the writing bigwig file
  covFile <- file.path(
    outdir,
    paste0(groupNamei, "-TileSize-",tileSize,"-normMethod-",normMethod,".bw")
  )
  
  covList <- lapply(X = seq_along(availableChr), FUN = function(k) {
    fragik <- fragi[seqnames(fragi) == availableChr[k],]
    tilesk <- tiles[BiocGenerics::which(S4Vectors::match(seqnames(tiles), availableChr[k], nomatch = 0) > 0)]
    if (length(x = fragik) == 0) {
      tilesk$reads <- 0
      # If fragments
    } else {
      # N Tiles
      nTiles <- chromLengths[availableChr[k]] / tileSize
      # Add one tile if there is extra bases
      if (nTiles%%1 != 0) {
        nTiles <- trunc(x = nTiles) + 1
      }
      # Create Sparse Matrix
      matchID <- S4Vectors::match(mcols(fragik)$name, cellGroupi)
      
      # For each tiles of this chromosome, create start tile and end tile row,
      # set the associated counts matching with the fragments
      mat <- sparseMatrix(
        i = c(trunc(x = start(x = fragik) / tileSize),
              trunc(x = end(x = fragik) / tileSize)) + 1,
        j = as.vector(x = c(matchID, matchID)),
        x = rep(1, 2*length(x = fragik)),
        dims = c(nTiles, length(x = cellGroupi))
      )
      
      # Max count for a cells in a tile is set to cutoff
      if (!is.null(x = cutoff)){
        mat@x[mat@x > cutoff] <- cutoff
      }
      # Sums the cells
      mat <- Matrix::rowSums(x = mat)
      tilesk$reads <- mat
      # Normalization
      if (!is.null(x = normMethod)) {
        if (normMethod == "rc") {
          tilesk$reads <- tilesk$reads * 10^4 / length(fragi$name)
        } else if (normMethod == "ncells") {
          tilesk$reads <- tilesk$reads / length(cellGroupi)
        } else if (normMethod == "none") {
        } else {
          if (!is.null(x = normBy)){
            tilesk$reads <- tilesk$reads * 10^4 / sum(normBy[cellGroupi, 1])
          }
        }
      }
    }
    tilesk <- coverage(tilesk, weight = tilesk$reads)[[availableChr[k]]]
    tilesk
  })
  
  names(covList) <- availableChr
  covList <- as(object = covList, Class = "RleList")
  rtracklayer::export.bw(object = covList, con = covFile)
  return(covFile)
}

# Set a default value if an object is null
#
# @param x An object to set if it's null
# @param y The value to provide if x is null
# @return Returns y if x is null, otherwise returns x.
SetIfNull <- function(x, y) {
  if (is.null(x = x)) {
    return(y)
  } else {
    return(x)
  }
}
