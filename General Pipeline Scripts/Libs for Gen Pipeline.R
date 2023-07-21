#' Libs for Gen Pipeline
#' 6/19/23
#'
#' Written by: Benji Bromberg
#'
#' Purpose:
#'  The purpose of this script is to load all libraries that are used in the
#'  pipeline.
#'
#'  Note: Some libraries have conflicting function names and therefore the
#'  library name may need to be specified for proper function use. The most
#'  frequent case of conflict is when calling the dplyr::select() function.
#'
#'  Outdated packages/includes:
#'  source("api.R")

library(tidyverse)
library(dplyr)
library(readxl)
library(ggrepel)
library(pheatmap)
library(ashr)
library(stringr)
library(tidyr)
library(tibble)
library(utils)
library(R.utils)
library(DBI)
library(RSQLite)
library(httpgd)
library(lintr)

# Bioconductor packages
library(enrichplot)
library(GO.db)
library(DESeq2)
library(clusterProfiler)
library(edgeR)
library(vsn)
library(biomaRt)
library(DOSE)
library(AnnotationForge)
library(AnnotationDbi)
library(pathview)

# Attempting to load org.Sjaponicus.eg.db
tryCatch(
  {
    library(org.Sjaponicus.eg.db) # Bioconductor package
  },
  error = function(e) {
    stop("Error: There is no package called 'org.Sjaponicus.eg.db'. Please set
      working directory to /Mack RNA Seq/Japonicus/SJaponicusFromEggnog/ and
      install 'org.Sjaponicus.eg.db' using the command below.

      setwd(paste0(script_dir, '/Japonicus/SJaponicusFromEggnog/'))
      install.packages('org.Sjaponicus.eg.db/', repos = NULL, type = 'source')")
  }
)

# Attempting to load org.Tmaccoyyi.eg.db
tryCatch(
  {
    library(org.Tmaccoyii.eg.db) # Bioconductor package
  },
  error = function(e) {
    stop("Error: There is no package called 'org.Tmaccoyyi.eg.db'. Please set
      working directory to /Mack RNA Seq/Maccoyii/TMaccoyiiFromEggnog/ and
      install 'org.Tmaccoyyi.eg.db' using the command below.

      setwd(paste0(script_dir, '/Maccoyii/TMaccoyiiFromEggnog/'))
      install.packages('org.Tmaccoyii.eg.db/', repos = NULL, type = 'source')")
  }
)
