# Package R dependencies

dependencies <- c("dada2","phyloseq","Biostrings","ComplexHeatmap","ggtext","ggplot2","dplyr","pheatmap","reshape2",
                  "ape","microbiome","DT","vegan","ggthemes","ggsignif","GUniFrac", "tibble","tidyverse","magrittr",
                  "readxl","multcomp","mvtnorm", "ggpicrust2","phangorn","biomformat", "ggprism", "viridisLite", "viridis",
                  "survival","TH.data","MASS","ggraph","corncob","adegenet","ade4","KEGGREST","edgeR","readr")

# devtools 

if( !is.element("devtools",rownames(installed.packages() ) ) ){
  install.packages("devtools")
  install.packages("BiocManager")
}
library(devtools)

# Install missing packages

missingPackages <- function(pkg){
  if( !is.element(pkg,rownames(installed.packages() ) ) ){
    message(pkg, "-----> Package is not installed ")
    BiocManager::install(pkg)
  }
}

for(i in dependencies){
  missingPackages(i)
  library(i, character.only = TRUE)
}

if( !is.element("microViz",rownames(installed.packages() ) ) ){
  install.packages(
    "microViz",
    repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
  )
}
library(microViz)
