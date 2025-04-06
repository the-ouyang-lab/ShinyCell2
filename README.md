# ShinyCell2
`ShinyCell2` is an enhanced R package for creating interactive, lightweight, 
and shareable web applications to explore single-cell multi-omics and spatial 
transcriptomics data. Building on the simplicity of the original ShinyCell, 
ShinyCell2 introduces powerful new features tailored for modern single-cell 
modalities, including support for CITE-seq, scATAC-seq, and spatial 
transcriptomics. It integrates seamlessly with popular analysis tools like 
Seurat, Scanpy, Signac, and ArchR, and offers advanced visualizations such as 
zoom-enabled UMAPs, IGV-style peak tracks, and spatial plots—all easily 
customizable and deployable with minimal dependencies. Designed for both 
computational and experimental researchers, ShinyCell2 empowers intuitive 
exploration, cross-modality comparison, and statistical analysis of 
high-dimensional data without requiring extensive coding or setup. An example 
web app showcasing the visualization of spatial transcriptomics, scATAC-seq, 
and CITE-seq datasets can be accessed at https://shinycell.ouyanglab.com/

If you are using `ShinyCell2`, please cite the biorxiv preprint. The manuscript 
is currently under review.

Key features of `ShinyCell2` include:

- Supports seamless integration and visualization of single-cell multi-omics, 
peak-based, and spatial data for in-depth exploration.

- Converts data from leading analysis tools—including Scanpy, Seurat, Signac, 
and ArchR—into a lightweight, portable format suitable for local or web-based 
deployment, enabling open and accessible data sharing.

- Offers a rich collection of interactive, customizable, and publication-ready 
plots.

- Easily extendable, with support for user-defined visualizations through R for 
tailored analysis workflows.



# Table of Contents and Additional Information / Tutorials
This readme is broken down into the following sections:

- [Installation](#installation) on how to install `ShinyCell2`

- [Quick Start Guide](#quick-start-guide) to rapidly deploy a shiny app with 
  a few lines of code

- [Frequently Asked Questions](#frequently-asked-questions)

There are also additional information / tutorials as follows:

- [Additional information on new visualisations tailored for spatial / scATAC-seq / multiomics](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Additional information on enhanced visualisation features](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell2 app for scRNA-seq or CITE-seq data](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell2 app for spatial data](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell2 app for ArchR-based scATAC-seq data](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell2 app for Signac-based scATAC-seq data](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Tutorial for creating a ShinyCell app containing several sdatasets](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/2multi.html)

- [Tutorial for customising ShinyCell aesthetics](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)

- [Instructions on how to deploy ShinyCell apps online](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html)



# Installation
First, users can run the following code to check if the packages required by 
`ShinyCell2` exist and install them if required:
``` r
reqPkg = c("data.table", "Matrix", "hdf5r", "reticulate", "R.utils", 
           "ggplot2", "gridExtra", "glue", "readr", "future", "RColorBrewer")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}

# If you are using Seurat object as input (for scRNA, multiomics, spatial), 
#   you can install Seurat as follows:
# install.packages("Seurat")

# If you are using Signac object for scATAC-seq, you can install Signac as follows:
# install.packages("Signac")

# If you are using ArchR object as input for scATAC-seq, visit
#   https://github.com/GreenleafLab/ArchR for ArchR's installation instruction

# If you are using h5ad file as input, run code below
# reticulate::py_install("anndata")
```

Furthermore, on the system where the Shiny app will be deployed, users can run 
the following code to check if the packages required by the Shiny app exist 
and install them if required:
``` r
reqPkg = c("shiny", "shinyhelper", "data.table", "Matrix", "DT", "magrittr", 
           "ggplot2", "ggrepel", "hdf5r", "ggdendro", "gridExtra", "ggpubr")
newPkg = reqPkg[!(reqPkg %in% installed.packages()[,"Package"])]
if(length(newPkg)){install.packages(newPkg)}
```

If one is deploying scATAC-seq / peak-based data, bwtools have to be installed 
for the track plots via the command line:
``` r
# Install dependencies for trackplot (on command line)
git clone https://github.com/CRG-Barcelona/libbeato.git
cd ./libbeato
git checkout 0c30432
./configure
make
make install
    
cd ..
git clone https://github.com/CRG-Barcelona/bwtool.git
cd ./bwtool
./configure
make
make check
make install
```

`ShinyCell2` can then be installed from GitHub as follows:
``` r
devtools::install_github("the-ouyang-lab/ShinyCell2")
```



# Quick Start Guide
In short, the `ShinyCell2` package takes in an input single-cell object and 
generates a ShinyCell2 config `scConf` containing labelling and colour palette 
information for the single-cell metadata. The ShinyCell config and single-cell 
object are then used to generate the files and code required for the shiny app. 

In this example, we will use single-cell CITE-seq data in the form of a Seurat 
object containing 162,000 PBMC cells measured with 228 antibodies, taken from 
https://satijalab.org/seurat/articles/multimodal_reference_mapping.html
The Seurat object can be [downloaded here](
https://zenodo.org/records/15162323/files/multimodal_pbmc.rds?download=1).

A shiny app can be readily generated using the following code:
 
``` r
library(Seurat)
library(ShinyCell)

seu <- readRDS("multimodal_pbmc.rds")
scConf <- createConfig(seu)
makeShinyFiles(seu,scConf, shiny.prefix="sc1", shiny.dir="shinyApp/")
makeShinyCodes(shiny.title = "PBMC multiomics", shiny.prefix="sc1",
               shiny.dir="shinyApp/")
```

The generated shiny app can then be found in the `shinyApp/` folder (which is 
the default output folder). To run the app locally, use RStudio to open either 
`server.R` or `ui.R` in the shiny app folder and click on "Run App" in the top 
right corner. The shiny app can also be deployed online via online platforms 
e.g. [shinyapps.io](https://www.shinyapps.io/) and Amazon Web Services (AWS) 
or be hosted via Shiny Server. For further details, refer to 
[Instructions on how to deploy ShinyCell apps online](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/4cloud.html).

More details on the various visualisations in the `ShinyCell2` can be found in
[Additional information on new visualisations tailored for spatial / scATAC-seq / multiomics](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)
and [Additional information on enhanced visualisation features](
https://htmlpreview.github.io/?https://github.com/SGDDNB/ShinyCell/blob/master/docs/1aesthetics.html)



# Frequently Asked Questions
- Q: How much memory / storage space does `ShinyCell2` and the app consume?
  - A: The `ShinyCell2` app consumes very little memory and is meant to be a 
       heavy-duty app where multiple users can access the app simultaneously. 
       Unlike typical R objects, the entire gene expression matrix is stored 
       on disk and *not on memory* via the hdf5 file system. Also, the hdf5 
       file system offers superior file compression and takes up less storage 
       space than native R file formats such as rds / Rdata files.
  - A: It should be noted that a large amount of memory is required when 
       *building* the `ShinyCell2` app. This is because the whole single-cell 
       object has to be loaded onto memory and additional memory is required to 
       generate the required files. From experience, a typical laptop with 8GB 
       RAM can handle datasets around 30k cells while 16GB RAM machines can 
       handle around 60k-70k cells for scRNA-seq data. More memory will be 
       required for multiomics and especially scATAC-seq data. As a rule of 
       thumb, if you are able to perform analysis e.g. dimension reduction on 
       the machine, the machine should be able to build the `ShinyCell2` app.
       
- Q: I have both RNA and integrated data in my Seurat object. How do I specify 
which gene expression assay to plot in the Shiny app?
  - A: Unlike the original ShinyCell, `ShinyCell2` now supports multiple assays 
       within a single `ShinyCell2` app. Thus, both the RNA and integrated data 
       will be incorporated and users can choose to visualise either assays or 
       even compare their expression in the `ShinyCell2` app.


