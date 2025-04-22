#' Generate code files required for shiny app
#'
#' Generate code files required for shiny app for both single-dataset and 
#' multi-dataset scenarios. Specifically, two R scripts will be generated, 
#' namely \code{server.R} and \code{ui.R}. Note that \code{makeShinyFiles} has 
#' to be ran prior to make the necessary data files for each dataset included. 
#' The prefix used in \code{makeShinyFiles} have to be then supplied in this 
#' function.
#'
#' @param shiny.title specify the overall title for shiny app
#' @param shiny.footnotes text for shiny app footnote. When given as a list, 
#'   citation can be inserted by specifying author, title, journal, volume, 
#'   page, year, doi, link. See example below. 
#' @param shiny.prefix specify file prefix for each dataset. Must match the 
#'   prefix used in \code{makeShinyFiles}
#' @param shiny.headers specify the tab header names for each dataset. Length 
#'   must match that of \code{shiny.prefix}. Note that this is ignored if 
#'   there is only one dataset
#' @param shiny.dir specify directory to create the shiny app in
#' @param defPtSiz specify default point size for single cells. For example, a 
#'   smaller size can be used if you have many cells in your dataset. A single 
#'   value can be specified to set the point size for all datasets. Otherwise,
#'   users have to specify one value for each dataset
#' @param ganalytics Google analytics tracking ID (e.g. "UA-123456789-0")
#'
#' @return server.R and ui.R required for shiny app
#'
#' @author John F. Ouyang
#'
#' @import data.table readr
#'
#' @examples
#' # Example citation
#' citation = list(
#'   author  = "Liu X., Ouyang J.F., Rossello F.J. et al.",
#'   title   = "",
#'   journal = "Nature",
#'   volume  = "586",
#'   page    = "101-107",
#'   year    = "2020", 
#'   doi     = "10.1038/s41586-020-2734-6",
#'   link    = "https://www.nature.com/articles/s41586-020-2734-6")
#' makeShinyCodes(shiny.title = "scRNA-seq shiny app", shiny.footnotes = "",
#'                shiny.prefix = c("sc1", "sc2"), defPtSiz = c(1.25, 1.5),
#'                shiny.headers = c("dataset1", "dataset2"),
#'                shiny.dir = "shinyApp/")
#'                
#' @export
makeShinyCodes <- function(shiny.title, shiny.footnotes = "",
                           shiny.prefix, shiny.headers, shiny.dir, 
                           defPtSiz = 1.25, ganalytics = NA){
  ### Checks
  if(length(shiny.prefix) > 1){
    if(length(shiny.prefix) != length(shiny.headers)){
      stop("length of shiny.prefix and shiny.headers does not match!")
    }
  }
  if(length(shiny.prefix) != length(defPtSiz)){
    defPtSiz = rep(defPtSiz[1], length(shiny.prefix))
  }
  defPtSiz2 = as.character(defPtSiz*2)
  defPtSiz = as.character(defPtSiz)
  
  ### Check if files exist
  for(i in shiny.prefix){
    if(!file.exists(paste0(shiny.dir, "/", i, "conf.rds"))){
      stop(paste0("files missing for ", i, " dataset!"))
    }
  }
  
  ### Write code for shinyFunc.R
  fname = paste0(shiny.dir, "/shinyFunc.R")
  readr::write_file(wrShFunc(), file = fname)
  for(i in shiny.prefix){
    if(file.exists(paste0(shiny.dir, "/", i, "bw.rds"))){
      readr::write_file(wrShFuncT1(), file = fname, append = TRUE)
      srcPath <- system.file("extdata", "trackplot.R", package = "ShinyCell2")
      tarPath <- paste0(shiny.dir,"/")
      file.copy(srcPath, tarPath)
    }
  }

  ### Write code for server.R
  # Function will auto detect if ATAC or spatial is present
  fname = paste0(shiny.dir, "/server.R")
  readr::write_file(wrSVlib(), file = fname)
  for(i in shiny.prefix){
    readr::write_file(wrSVload(i), file = fname, append = TRUE)
    if(file.exists(paste0(shiny.dir, "/", i, "image.rds"))){
      readr::write_file(wrSVloadS1(i), file = fname, append = TRUE)}
    if(file.exists(paste0(shiny.dir, "/", i, "bw.rds"))){ 
      readr::write_file(wrSVloadT1(i), file = fname, append = TRUE)}
  }
  readr::write_file(wrSVpre(), file = fname, append = TRUE)
  for(i in shiny.prefix){
    if(file.exists(paste0(shiny.dir, "/", i, "image.rds"))){
      readr::write_file(wrSVmainS1(i), file = fname, append = TRUE)
      readr::write_file(wrSVmainS2(i), file = fname, append = TRUE)}
    if(file.exists(paste0(shiny.dir, "/", i, "bw.rds"))){ 
      readr::write_file(wrSVmainT1(i), file = fname, append = TRUE)}
    readr::write_file(wrSVmainA1(i), file = fname, append = TRUE)
    readr::write_file(wrSVmainA2(i), file = fname, append = TRUE)
    readr::write_file(wrSVmainA3(i), file = fname, append = TRUE)
    readr::write_file(wrSVmainB1(i), file = fname, append = TRUE)
    readr::write_file(wrSVmainB2(i), file = fname, append = TRUE)
    readr::write_file(wrSVmainB3(i), file = fname, append = TRUE)
  }
  readr::write_file(wrSVpost(), file = fname, append = TRUE)

  ### Write code for ui.R
  fname = paste0(shiny.dir, "/ui.R")
  readr::write_file(wrUIlib(), file = fname)
  for(i in shiny.prefix){
    readr::write_file(wrUIload(i), file = fname, append = TRUE)
    if(file.exists(paste0(shiny.dir, "/", i, "bw.rds"))){ 
      readr::write_file(wrUIloadT1(i), file = fname, append = TRUE)}
  }
  readr::write_file(wrUIpre(shiny.title, ganalytics), file = fname, append = TRUE)
  if(length(shiny.prefix) == 1){
    if(file.exists(paste0(shiny.dir, "/", shiny.prefix, "image.rds"))){
      readr::write_file(wrUImainS1(shiny.prefix, defPtSiz2[1]), file = fname, append = TRUE)
      readr::write_file(wrUImainS2(shiny.prefix, defPtSiz2[1]), file = fname, append = TRUE)}
    if(file.exists(paste0(shiny.dir, "/", shiny.prefix, "bw.rds"))){
      readr::write_file(wrUImainT1(shiny.prefix), file = fname, append = TRUE)}
    readr::write_file(wrUImainA1(shiny.prefix, defPtSiz[1]), file = fname, append = TRUE)
    readr::write_file(wrUImainA2(shiny.prefix, defPtSiz[1]), file = fname, append = TRUE)
    readr::write_file(wrUImainA3(shiny.prefix, defPtSiz[1]), file = fname, append = TRUE)
    readr::write_file(wrUImainB1(shiny.prefix, defPtSiz[1]), file = fname, append = TRUE)
    readr::write_file(wrUImainB2(shiny.prefix), file = fname, append = TRUE)
    readr::write_file(wrUImainB3(shiny.prefix), file = fname, append = TRUE)
    readr::write_file(glue::glue(', \n'), append = TRUE, file = fname)
  } else {
    for(i in seq_along(shiny.prefix)){
      hhh = shiny.headers[i]
      readr::write_file(glue::glue('navbarMenu("{hhh}",'), file = fname, append = TRUE)
      if(file.exists(paste0(shiny.dir, "/", shiny.prefix[i], "image.rds"))){
        readr::write_file(wrUImainS1(shiny.prefix[i], defPtSiz2[i]), file = fname, append = TRUE)
        readr::write_file(wrUImainS2(shiny.prefix[i], defPtSiz2[i]), file = fname, append = TRUE)}
      if(file.exists(paste0(shiny.dir, "/", shiny.prefix[i], "bw.rds"))){ 
        readr::write_file(wrUImainT1(shiny.prefix[i]), file = fname, append = TRUE)}
      readr::write_file(wrUImainA1(shiny.prefix[i], defPtSiz[i]), file = fname, append = TRUE)
      readr::write_file(wrUImainA2(shiny.prefix[i], defPtSiz[i]), file = fname, append = TRUE)
      readr::write_file(wrUImainA3(shiny.prefix[i], defPtSiz[i]), file = fname, append = TRUE)
      readr::write_file(wrUImainB1(shiny.prefix[i], defPtSiz[i]), file = fname, append = TRUE)
      readr::write_file(wrUImainB2(shiny.prefix[i]), file = fname, append = TRUE)
      readr::write_file(wrUImainB3(shiny.prefix[i]), file = fname, append = TRUE)
      readr::write_file(glue::glue('), \n\n\n'), append = TRUE, file = fname)
    }
  }
  readr::write_file(wrUIpost(shiny.footnotes), file = fname, append = TRUE)
  
  ### Write code for google-analytics.html
  if(!is.na(ganalytics)){
    fname = paste0(shiny.dir, "/google-analytics.html")
    readr::write_file(wrUIga(ganalytics), file = fname)
  }
  
}


