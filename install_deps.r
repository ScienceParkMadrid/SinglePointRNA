# These libraries/packages (Ubuntu) are needed to build some of the dependencies
# git
# libpq-dev
# libhdf5-dev
# liblzma-dev
# libbz2-dev
# libglpk-dev
# libfftw3-3


# dependencies
packages <- c(
    "shinydashboard",
    "shinyFiles",
    "dplyr",
    "ggplot2",
    "grid",
    "gridExtra",
    "plotly",
    "ggalluvial",
    "openxlsx",
    "zoo",
    "readxl",
    "Seurat",
    "hdf5r",
    "HH",
    "stringr",
    "sctransform",
    "scCATCH",
    "BiocParallel",
    "glmGamPoi"
)

# Install BiocManager first
if(!require('BiocManager',character.only = TRUE)) install.packages("BiocManager", dependencies = TRUE)

# Install the dependencies
lapply(packages, function(x) if(!require(x,character.only = TRUE)) BiocManager::install(x, dependencies = TRUE, ask = FALSE))
