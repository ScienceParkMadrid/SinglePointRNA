# SinglePoint RNA

SinglePoint RNA is a shiny-based R application that provides a graphic interface for different publicly available tools to analyze single cell RNA-seq data. 

The aim of SinglePoint RNA is to provide an accessible and transparent tool set to researchers that allows them to perform detailed and custom analysis of their data autonomously. 
SinglePointRNA is structured in a context-driven framework that prioritizes providing the user with solid qualitative guidance at each step of the analysis process and interpretation of the results.
Additionally, the rich user guides accompanying the software are intended to serve as a point of entry for users to learn more about computational techniques applied to single cell data analysis.


## Installation

### Requirements 

* [R](https://cran.r-project.org/)
* [RStudio](https://posit.co/download/rstudio-desktop/)

SinglePoint RNA depends on the following R packages:

* [shiny](https://shiny.rstudio.com/)
* [shinydashboard](https://cran.r-project.org/package=shinydashboard)
* [shinyFiles](https://cran.r-project.org/package=shinyFiles)
* [dplyr](https://cran.r-project.org/package=dplyr)
* [ggplot2](https://cran.r-project.org/package=ggplot2) 
* [grid](https://cran.r-project.org/package=grid) 
* [gridExtra](https://cran.r-project.org/package=gridExtra) 
* [plotly](https://cran.r-project.org/package=plotly)
* [ggalluvial](https://cran.r-project.org/package=ggalluvial)
* [openxlsx](https://cran.r-project.org/package=openxlsx)
* [zoo](https://cran.r-project.org/package=zoo)
* [readxl](https://cran.r-project.org/package=readxl)
* [Seurat](https://cran.r-project.org/package=Seurat)
* [hdf5r](https://cran.r-project.org/package=hdf5r)
* [HH](https://cran.r-project.org/package=HH)
* [stringr](https://cran.r-project.org/package=stringr)
* [sctransform](https://cran.r-project.org/package=sctransform)
* [scCATCH](https://cran.r-project.org/package=scCATCH)

* [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html)
* [glmGamPoi](https://bioconductor.org/packages/release/bioc/html/glmGamPoi.html)

You can install the dependencies running the `install_deps.r` script:

```console
$ Rscript install_deps.r
```

Memory and processing requirements will vary depending on the size of the datasets to analyzed. As a reference, for a 5,000 cell dataset, a minimum of 8GB and 4 cores is recommended.

### Running SinglePoint RNA

After installing R, RStudio, and the R libraries stated above:

* Download the folder "SinglePointRNA" from this repository. 
* Open the file "SinglePointRNA/app.R" with Rstudio.
* Click the "Run App" button at the top of the source pane.


It's also possible to run SinglePoint RNA in a remote server using [Shiny Server](https://posit.co/products/open-source/shinyserver/)

## Citation

[SinglePointRNA, an user-friendly application implementing single cell RNA-seq analysis software](https://arxiv.org/abs/2305.00008)
