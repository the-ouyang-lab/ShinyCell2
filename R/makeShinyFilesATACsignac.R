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
#' @param Genome assembly used for annotation. Must be one of: hg19, hg38, mm10, mm9
#' @param chrom_assay_name Name of chromatin assay (Signac)
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
#' @importFrom IRanges IRanges start end
#' @importFrom GenomicRanges GRanges slidingWindows mcols seqnames
#' @importClassesFrom IRanges IRanges
#' @importClassesFrom GenomicRanges GRanges
#' 
#' @examples
#' makeShinyFilesATACsignac(sig, scConf, 
#'                          shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyFilesATACsignac <- function(
  obj, scConf, bigWigGroup, shiny.prefix = "sc1", shiny.dir = "shinyApp/", genome, chrom_assay_name){

  # Check where the seqinfo is stored and copy if necessary
  if (is.null(obj@assays[[chrom_assay_name]]@seqinfo)){
      obj@assays[[chrom_assay_name]]@seqinfo <- obj@assays[[chrom_assay_name]]@annotation@seqinfo
  }
    
  # Make individual bigwig files
  sc1conf = scConf
  for(iGrp in bigWigGroup){
    if(iGrp %in% sc1conf[grp == TRUE]$ID){
      if(!dir.exists(paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp))){
        dir.create(paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp))
      }
      ExportGroupBW(
        obj, assay = chrom_assay_name, group.by = iGrp,
        idents = NULL, normMethod = "RC",              # defaults
        tileSize = 100, minCells = 5, cutoff = NULL,   # defaults
        chromosome = GenomeInfoDb::standardChromosomes(obj),
        outdir = paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp), verbose = FALSE)
      system2(command = "rm", args = c(
        "-rf", paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp,"/*.bed")))
      bwFile <- list.files(path = paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp), pattern = "*.bw")
      for(iFile in bwFile){
        iFile2 <- iFile; if(grepl("^X[0-9]+", iFile)){iFile2 <- gsub("^X", "", iFile)}
        iFile2 <- gsub("-TileSize-100-normMethod-rc.bw$", ".bw", iFile2)
        iFile1 <- paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp,"/",iFile)
        iFile2 <- paste0(shiny.dir,"/",shiny.prefix,"bw_",iGrp,"/",iFile2)
        if(iFile1 != iFile2){system2(command = "mv", args = c(iFile1, iFile2))}
      }
    } else {
      warning(paste0(iGrp, " is not present in bigWigGroup. Skipping bigWig generation!"))
    }
  }
  # Copy over gtf
  scGenome <- genome
  # scGenome <- obj@assays[["peaks"]]@annotation@seqinfo@genome[1]
  # scGenome <- unlist(strsplit(getGenome(obj), "\\."))
  # scGenome <- scGenome[grepl("hg|mm", scGenome)]
  srcPath <- system.file("extdata", paste0(scGenome,".refGene.gtf.gz"), package = "ShinyCell2")
  tarPath <- paste0(shiny.dir,"/",shiny.prefix,"bw.gtf.gz")
  file.copy(srcPath, tarPath)

  # Generate geneIndex and chrom.sizes
  chrSize <- data.table(chr = as.character(seqnames(obj@assays[[chrom_assay_name]]@seqinfo)), 
                        size = as.numeric(seqlengths(obj@assays[[chrom_assay_name]]@seqinfo)))
  chrSize <- chrSize[chr %in% GenomeInfoDb::standardChromosomes(obj)]
  chrGene <- obj@assays[[chrom_assay_name]]@annotation
  chrGene <- data.table(gene = chrGene$gene_name,
                        chr = as.character(seqnames(chrGene)), 
                        srt = start(chrGene),
                        end = end(chrGene), type = chrGene$type)
  chrGene <- chrGene[type == "exon"][, .(chr = unique(chr),
                                         srt = min(srt),
                                         end = max(end)), by = "gene"]
  sc1bw <- list(chrSize = chrSize, chrGene = chrGene,
                genome = scGenome)
  
  ### Saving objects
  saveRDS(sc1bw,   file = paste0(shiny.dir, "/", shiny.prefix, "bw.rds"))
  return(paste0(shiny.prefix, " bw files generated!"))
}


