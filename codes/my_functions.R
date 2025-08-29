# Metadata function ------------------------------------------------------------
# function to get metadata from GSE series matrix.

my_metadata_function <- function(my_GSE){
  
  gset <- getGEO(my_GSE, GSEMatrix =TRUE, getGPL=FALSE)
  gset <- gset[[1]]
  
  suppl_names = getGEOSuppFiles(my_GSE, makeDirectory = FALSE, baseDir = getwd(),
                                fetch_files = FALSE, filter_regex = NULL)

  sampleInfo = pData(gset)
  if ("disease state:ch1" %in% colnames(sampleInfo)) {
    #select cols
    sampleInfo2 = dplyr::select(sampleInfo, 
                         "title",
                         "age:ch1", 
                         "disease state:ch1", 
                         "gender:ch1", 
                         "geo_accession",
                         "supplementary_file")
    
    # rename cols
    sampleInfo2 = dplyr::rename(sampleInfo2,
                         Group = "disease state:ch1",
                         Age = "age:ch1",
                         Sex = "gender:ch1",
                         GSM = "geo_accession")
    
  } else {
    #select cols
    sampleInfo2 = dplyr::select(sampleInfo, 
                         "title",
                         "age:ch1", 
                         "source_name_ch1", 
                         "gender:ch1", 
                         "geo_accession",
                         "supplementary_file")
    
    # rename cols
    sampleInfo2 = dplyr::rename(sampleInfo2,
                         Group = "source_name_ch1",
                         Age = "age:ch1",
                         Sex = "gender:ch1",
                         GSM = "geo_accession") 
  }
  
  sampleInfo2$Group = str_to_lower(sampleInfo2$Group)
  sampleInfo2 = dplyr::mutate(sampleInfo2, .keep = "unused",
                              Group = gsub(x = sampleInfo2$Group, 
                                           pattern = " serum.*| sample.*", 
                                           replacement = ""))
  sampleInfo2 = dplyr::mutate(sampleInfo2, Group_original = Group)
  
  controlsindex = (grep(pattern =  "control", x = sampleInfo2$Group))
  sampleInfo2$Group[controlsindex] = "control"
  MSindex = (grep(pattern =  "sclerosis", x = sampleInfo2$Group))
  sampleInfo2$Group[MSindex] = "multiple sclerosis"
  PDindex = (grep(pattern =  "parkinson", x = sampleInfo2$Group))
  sampleInfo2$Group[PDindex] = "parkinson's disease"  
  ADindex = (grep(pattern =  "cognitive", x = sampleInfo2$Group))
  sampleInfo2$Group[ADindex] = "alzheimer's disease"
  
  sampleInfo2$Sex = str_to_lower(sampleInfo2$Sex)
  sampleInfo2$Age = str_to_lower(sampleInfo2$Age)
  #  sampleInfo2 = sampleInfo2 %>% dplyr::mutate(
  #  Sample_ID = str_extract(string = sampleInfo2$title, pattern = "[:digit:]*$"))
  sampleInfo2 = dplyr::mutate(sampleInfo2,
                              filename = gsub(x = sampleInfo2$supplementary_file, 
                                              pattern = ".*GSM[[:digit:]_]{7,8}|[.]gpr[.]gz$", 
                                              replacement = ""))
  sampleInfo2 = sampleInfo2 %>% dplyr::mutate(
    PA_ID = str_extract(string = sampleInfo2$filename, pattern = "PA[:digit:]{4,}"))
  
  #sampleInfo2 = sampleInfo2 %>% dplyr::arrange(Group, Sex, Age, PA_ID)
  rownames(sampleInfo2) <- NULL
  sampleInfo2$GSE = my_GSE
  
  sampleInfo2 = sampleInfo2[,c("title", "filename", "PA_ID", 
                               "Age", "Sex", "Group", "GSE", "GSM", 
                               "Group_original", "supplementary_file")]
  
  #nome = paste0("./data/Nagele_sampleinfo/", today(),"_" , my_GSE, ".csv")
  #write_excel_csv(x = sampleInfo2, file = nome)
  
  return(sampleInfo2)
}


# paweR aggregateElist ---------------------------------------------------------
#' Aggregating elist intensities so that only unique proteins remain
#'
#' @param elist An object with original intensity values.
#' @param by A name of the column, by which values would be aggregated.
#' @param method A name of the method for aggregation. Possible options are:
#' "mean", "median", "min", "max" and "none".
#'
#' @return targets object of class data.frame containing ArrayID and FileNames
#'
aggregateElist <- function(elist, by, method) {
  elistAggregated <- elist
  
  if (method != 'none') {
    elistAggregated$E <-
      stats::aggregate(elistAggregated$E, list(elistAggregated$genes[,by]), method)
    rownames(elistAggregated$E) <- elistAggregated$E$Group.1
    elistAggregated$E$Group.1 <- NULL
    
    # Aggregate annotation about proteins
    genes <- elistAggregated$genes
    splitGenes <- split(genes, genes[by])
    
    genes <- do.call(rbind, lapply(splitGenes, function(x) {
      x$Column[1] <- paste(x$Column, collapse = ';')
      #x$ID[1] <-
      #  paste(sub('C.*', 'C', x$ID[1]), paste(sub('.*C', '', x$ID), collapse = ';'), sep = '')
      return(x[1,])
    }))
    
    # TBA: Catch an error if they are different
    if (nrow(genes) == nrow(elistAggregated$E)) {
      elistAggregated$genes <- genes
    } else {
      # FOR DEBUG
      #message(nrow(genes))
      #message(nrow(elistAggregated$E))
      stop("Error: aggregated dimensions of matrix E and genes do not match. Interrupt.")
    }
  }
  return(elistAggregated)
}

