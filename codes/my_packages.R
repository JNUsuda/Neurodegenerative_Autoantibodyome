BiocManager::install(c("GEOquery", "limma", "sva", "DOSE", "clusterProfiler",
                       "rrvgo", "org.Hs.eg.db", "ComplexHeatmap"), force = F)

suppressMessages({
  
  load_lib <- c("tidyverse", "highcharter", "BiocManager", "forcats", "stringr",
                "ggrepel", "readr", "survminer", "pheatmap", "readxl", "svglite", 
                "ggridges", "ggvenn",  "corrplot", "factoextra", "berryFunctions", 
                "network", "GGally", "scales", "sna", "RColorBrewer", #"ggnet", 
                "ggpattern",
                "GEOquery", "limma",  "sva","DOSE","clusterProfiler", #bioconductor
                 "rrvgo", "org.Hs.eg.db", "ComplexHeatmap" ) #bioconductor

  install_lib <- load_lib[!(load_lib %in% installed.packages())] # check package
  if (length(install_lib)) for (i in install_lib) install.packages(i) # install
  
  cat("Loaded Packages:\n")
  print(sapply(load_lib, require, character = TRUE)) # load
  cat("\n\n")
  
})
