#' Generate data files required for shiny app (ArchR object)
#'
#' Generate data files required for shiny app, specifically scATAC-seq data
#' Six files will be generated, namely 
#' (i) the shinycell config \code{prefix_conf.rds}, 
#' (ii) the single-cell metadata \code{prefix_meta.rds}, 
#' (iii) the single-cell assays \code{prefix_assay_X.h5}, 
#' (iv) the feature mapping object config \code{prefix_gene.rds}, 
#' (v) the dimension reduction embeddings \code{prefix_dimr.rds} and 
#' (vi) the defaults for the Shiny app \code{prefix_def.rds} and
#' (vii) the bigwig files for trackplot \code{prefix_bw_GRP}. 
#' A prefix is specified for each set of files to allow for multiple 
#' single-cell datasets in a single Shiny app.
#'
#' @param obj input ArhcR object
#' @param scConf shinycell config data.table
#' @param bigWigGroup categorical group in ArhcR meta.data to group cells by for 
#'   the generation of bigWig files for track plot. Default is NA which does 
#'   not generate any bigWig files.
#' @param assay assay(s) in ArhcR object to use. Multiple assays can now be 
#'   incorporated and all assays are used by default (with the first assay 
#'   being the default assay), which must match one of the following:
#'   \itemize{
#'     \item{ArchR objects}: "TileMatrix" or "GeneScoreMatrix" or 
#'       "GeneIntegrationMatrix" or "PeakMatrix" or "MotifMatrix",
#'       default is "PeakMatrix"
#'   }
#' @param dimred.to.use specify the dimension reduction to use. Default is to 
#'   use all except LSI 
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#' @param default.gene1 specify primary default feature (peak or gene or TF) 
#'   to show, which must be present in the default assay
#' @param default.gene2 specify secondary default feature (peak or gene or TF) 
#'   to show, which must be present in the default assay
#' @param default.multigene character vector specifying default features to 
#'   show in bubbleplot / heatmap, which be present in the default assay
#' @param default.dimred character vector specifying the two default dimension 
#'   reductions. Default is to use UMAP if not TSNE embeddings
#' @param ... extra arguments to supply to ArchR::getGroupBW
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate
#'
#' @examples
#' makeShinyFilesATACarchr(ArchR, scConf, 
#'                         shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyFilesATACarchr <- function(
  obj, scConf, bigWigGroup = NA, assay = NA, dimred.to.use = NA,
  shiny.prefix = "sc1", shiny.dir = "shinyApp/",
  default.gene1 = NA, default.gene2 = NA, default.multigene = NA, 
  default.dimred = NA, ...){
  ### Preprocessing and checks
  if(class(obj) != "ArchRProject"){stop("obj is not of ArchR format!")}
  # Generate defaults for assay
  if(is.na(assay[1])){
    assay = ArchR::getAvailableMatrices(ArchRProj = obj)
    assay = setdiff(assay, c("TileMatrix","PeakMatrix"))
    if(length(assay) == 0){stop("No assays are present! Unable to proceed!")}
    if("GeneScoreMatrix" %in% assay){
      assay = c("GeneScoreMatrix", setdiff(assay, "GeneScoreMatrix"))
    }
  }
  # defgene = Seurat::VariableFeatures(obj)[1:10] # PROBLEM cannot define defgene
  if(is.na(dimred.to.use[1])){
    dimred.to.use = names(obj@embeddings@listData)
  }
  sc1meta = data.table(cellID = obj$cellNames, as.data.frame(obj@cellColData))
    
  # Check default.gene1 / default.gene2 / default.multigene
  # if(is.na(default.gene1[[1]])){default.gene1 = defGenes[1]}
  # if(is.na(default.gene2[[1]])){default.gene2 = defGenes[2]}
  # if(is.na(default.multigene[[1]])){default.multigene = defGenes}
  # default.dimred = default.dimred[1]
  # if(is.na(default.dimred)){default.dimred = dimred.to.use[1]}
  

  
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
  for(iDR in dimred.to.use){
    drMat = as.matrix(ArchR::getEmbedding(ArchRProj = obj, embedding = iDR, returnDF = TRUE))
    colnames(drMat) <- paste0(iDR, c("1", "2"))
    if(ncol(drMat) > 2){drMat = drMat[, 1:2]}  # Take first two comps only
    drMat = drMat[sc1meta$cellID, ]            # Ensure ordering
    sc1dimr[[iDR]] = drMat
  }
  
  # Make prefix_assay_X.h5 / prefix_gene.rds
  sc1gene = list()
  if(!dir.exists(shiny.dir)){dir.create(shiny.dir)}
  for(iAssay in assay){
    filename = paste0(shiny.dir, "/", shiny.prefix, "assay_", iAssay, ".h5")
    tmpG = makeH5fromArchR(obj, sc1meta, filename = filename, assay = iAssay)
    tmpOut = seq_along(tmpG); names(tmpOut) = tmpG
    tmpOut = tmpOut[order(names(tmpOut))]
    tmpOut = tmpOut[order(nchar(names(tmpOut)))]
    sc1gene[[iAssay]] = tmpOut
  }
  
  # Make bigwig files
  if(!is.na(bigWigGroup[1])){
    # Make individual bigwig files
    for(iGrp in bigWigGroup){
      if(iGrp %in% sc1conf[grp == TRUE]$ID){
        if(!dir.exists(paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp))){
          dir.create(paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp))
        }
        bwPath <- ArchR::getGroupBW(ArchRProj = obj, groupBy = iGrp, 
                                    normMethod = "ReadsInTSS", tileSize = 100, ...)
        system2(command = "mv", args = c(bwPath, paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp)))
        bwFile <- list.files(path = paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp), pattern = "*.bw")
        for(iFile in bwFile){
          iFile2 <- iFile; if(grepl("^X[0-9]+", iFile)){iFile2 <- gsub("^X", "", iFile)}
          iFile2 <- gsub("-TileSize-100-normMethod-ReadsInTSS-ArchR.bw$", ".bw", iFile2)
          iFile1 <- paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp,"/",iFile)
          iFile2 <- paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp,"/",iFile2)
          if(iFile1 != iFile2){system2(command = "mv", args = c(iFile1, iFile2))}
        }
      } else {
        warning(paste0(iGrp, " is not present in bigWigGroup. Skipping bigWig generation!"))
      }
    }
    # Copy over gtf
    scGenome <- unlist(strsplit(getGenome(obj), "\\."))
    scGenome <- scGenome[grepl("hg|mm", scGenome)]
    srcPath <- system.file("extdata", paste0(scGenome,".refGene.gtf.gz"), package = "ShinyCell2")
    tarPath <- paste0(shiny.dir,"/",shiny.prefix,"bw.gtf.gz")
    file.copy(srcPath, tarPath)
    # Generate geneIndex and chrom.sizes
    chrSize <- ArchR::getChromSizes(obj)
    chrSize <- data.table(chr = as.character(seqnames(chrSize)), 
                          size = end(chrSize))
    chrGene <- ArchR::getGenes(obj)
    chrGene <- data.table(gene = chrGene$symbol,
                          chr = as.character(seqnames(chrGene)), 
                          srt = start(chrGene),
                          end = end(chrGene))
    sc1bw <- list(chrSize = chrSize, chrGene = chrGene,
                  genome = scGenome)
    saveRDS(sc1bw,   file = paste0(shiny.dir, "/", shiny.prefix, "bw.rds"))
  }
  
  # Check (again) default.gene1 / default.gene2 / default.multigene
  defGene1 = list()
  defGene2 = list()
  defGeneM = list()
  for(iAssay in assay){
    tmp0 = names(sc1gene[[iAssay]][which(sc1gene[[iAssay]] %in% 1:10)])
    # default.gene1
    if(is.list(default.gene1) & (length(default.gene1) >= length(assay))){
      tmp = default.gene1[[iAssay == assay]][1]
    } else {
      tmp = default.gene1[[1]][1]
    }
    if(tmp %in% names(sc1gene[[iAssay]])){
      defGene1[[iAssay]] = tmp
    } else {
      defGene1[[iAssay]] = tmp0[1]
    }
    # default.gene2
    if(is.list(default.gene2) & (length(default.gene2) >= length(assay))){
      tmp = default.gene2[[iAssay == assay]][1]
    } else {
      tmp = default.gene2[[1]][1]
    }
    if(tmp %in% names(sc1gene[[iAssay]])){
      defGene2[[iAssay]] = tmp
    } else {
      defGene2[[iAssay]] = tmp0[2]
    }
    # default.multigene
    if(is.list(default.multigene) & (length(default.multigene) >= length(assay))){
      tmp = default.multigene[[iAssay == assay]]
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
  sc1def$assay = assay        # First one is default
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


