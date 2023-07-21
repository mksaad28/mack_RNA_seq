#' Make Org Packages From NCBI
#' 6/28/23
#' 
#' Written by: Benji Bromberg
#' 
#' Purpose:
#'  This script creates Org packages from using the makeOrgPackageNCBI function 
#'  found in the AnnotationForge package. Org packages are made for two species:
#'  Scomber japonicus and Scomber colias. 
#'  
#'  The script starts by loading all required packages and setting the working 
#'  directories for each species. For each species, the script calls the 
#'  makeOrgPackageFromNCBI function using the appropriate taxonomic information 
#'  for the species to create the package and then installs the package using 
#'  the install.packages function. The final output is two Org packages that can
#'  be used in further analysis.
#'  
#'  Notes: 
#'  - The makeOrgPackageFromNCBI function takes multiple hours to run. 
#'  - Only the final Org packages are uploaded to the GitHub due to file size 
#'  limits.

# Load all required packages
script_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd(script_dir)
source("General Pipeline Scripts/Libs for Gen Pipeline.R")

# Set timeout to be very long to prevent timeouts
options(timeout = 100000000000000)

# Make Scomber japonicus NCBI-based package ----
  # Set working directory
  script_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(paste0(script_dir, "/Japonicus/SJaponicusFromNCBI/"))

  # NOTE: Must uncomment to run (running takes multiple hours)
  
  # makeOrgPackageFromNCBI(version = "0.1",
  #                        author = 
  #                        "Benjamin Bromberg <bromberg.benji@gmail.com>",
  #                        maintainer = 
  #                        "Benjamin Bromberg <bromberg.benji@gmail.com>",
  #                        outputDir = ".",
  #                        tax_id = "13676",
  #                        genus = "Scomber",
  #                        species = "japonicus",
  #                        rebuildCache = TRUE)
  # install.packages("org.Sjaponicus.eg.db/", repos = NULL, type = "source")
  
# Make Scomber colias NCBI-based package ----
  # Set working directory
  script_dir <- dirname(dirname(rstudioapi::getActiveDocumentContext()$path))
  setwd(paste0(script_dir, "/Colias/SColiasFromNCBI/"))
  
  # NOTE: Must uncomment to run (running takes multiple hours)
  
  # makeOrgPackageFromNCBI(version = "0.1",
  #                        author =
  #                        "Benjamin Bromberg <bromberg.benji@gmail.com>",
  #                        maintainer =
  #                        "Benjamin Bromberg <bromberg.benji@gmail.com>",
  #                        outputDir = ".",
  #                        tax_id = "338315",
  #                        genus = "Scomber",
  #                        species = "colias",
  #                        rebuildCache = TRUE)
  # install.packages("org.Scolias.eg.db/", repos = NULL, type = "source")