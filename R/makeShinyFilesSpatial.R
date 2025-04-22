#' Generate data files required for shiny app (spatial data)
#'
#' Generate data files required for shiny app, specifically spatial data
#' A prefix is specified for each set of files to allow for multiple 
#' single-cell datasets in a single Shiny app.
#'
#' @param obj input Seurat (v3+) object or input file path for h5ad file
#' @param scConf shinycell config data.table
#' @param shiny.prefix specify file prefix 
#' @param shiny.dir specify directory to create the shiny app in
#'
#' @return data files required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table hdf5r reticulate hdf5r
#'
#' @examples
#' makeShinyFilesSpatial(seu, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/")
#'
#' @export
makeShinyFilesSpatial <- function(
  obj, scConf, shiny.prefix = "sc1", shiny.dir = "shinyApp/"){
  ### Preprocessing and checks
  if(!file.exists(paste0(shiny.dir, "/", shiny.prefix, "conf.rds"))){
    stop(paste0(shiny.prefix, "conf.rds file is missing! Have you ran makeShinyFilesGEX?"))
  }
  
  ### Start extraction
  sc1image = list()
  if(class(obj)[1] == "Seurat"){
    if(.hasSlot(obj, "images")){
      # coordinates
      sc1image$coord <- obj@images[[1]]@coordinates
      # uncropped background image
      bg_image <- GetImage(obj, mode = "raster") # background image
      # bg_grob <- rasterGrob(bg_image, width=unit(1,"npc"), height=unit(1,"npc"), 
      #                       interpolate = FALSE) # background image grob
      sc1image$bg_image <- bg_image
      # sc1image$bg_grob  <- bg_grob
      # for cropped image generation in the server.R script
      sc1image$lowres   <- obj@images[[1]]@scale.factors$lowres

    }
  } else {
    stop("Only Seurat objects are accepted!")
  }
  
  ### Saving objects
  saveRDS(sc1image, file = paste0(shiny.dir, "/", shiny.prefix, "image.rds"))
}


