multiOmicsViz <- function(sourceOmics, sourceOmicsName, chrome_sourceOmics,
                          targetOmicsList, targetOmicsName, chrome_targetOmics,
                          fdrThr, outputfile, nThreads=NULL, legend=TRUE) {
  
  outputfile <- paste(outputfile, ".pdf", sep="")
  
  if(class(sourceOmics) == "SummarizedExperiment") {
    sourceOmics <- assays(sourceOmics, n=1)
  } else {
    if(!is.matrix(sourceOmics) && !is.data.frame(sourceOmics)) {
      stop("Source omics data (e.g. CNA data) should be a R matrix, 
           data.frame or SummarizedExperiment object.")
    }
  }
  
  if(!is.character(sourceOmicsName)) {
    stop("sourceOmicsName is the name of the source omics data, which 
         should be a R character object.") 
  }
  
  if(!is.list(targetOmicsList)) {
    stop("Multiple target omics data (e.g. mRNA or protein data) 
         should be saved in the list object.")
  }
  
  if(length(targetOmicsList) > 5) {
    stop("targetOmicsList can only contain at most five omics data.")
  }
  
  for(i in seq_len(length(targetOmicsList))) {
    if(class(targetOmicsList[[i]]) == "SummarizedExperiment") {
      targetOmicsList[[i]] <- assays(targetOmicsList[[i]], n=1)
    } else {
      if(!is.matrix(targetOmicsList[[i]]) && !is.data.frame(targetOmicsList[[i]])) {
        stop("Each of all target omics data in the list (e.g. mRNA or 
             protein data) should be a R matrix, data.frame or 
             SummarizedExperiment object.")
      }
    }
  }
  
  if(!is.character(targetOmicsName)) {
    stop("targetOmicsName should be a R vector object 
         containing the name for each omics data in the targetOmicsList.")
  }
  
  if(length(targetOmicsList) != length(targetOmicsName)) {
    stop("targetOmicsName should have the same length 
         with targetOmicsList.")
  }   
  
  if(length(targetOmicsList) > 1 && is.null(nThreads)) {
    stop("Please input nThreads for the parallel computing.")
  }
  
  if(length(targetOmicsList) > 1 && !is.null(nThreads)) {
    if(nThreads > length(targetOmicsList)) {
      stop("nThreads should be at most the length of targetOmicsList.")
    }
  }
  
  ########find the intersect genes among all omics data######
  intG <- c()
  for(i in seq_len(length(targetOmicsList))) {
    if(i == 1) {
      intG <- rownames(targetOmicsList[[i]])
    } else {
      intG <- intersect(intG, rownames(targetOmicsList[[i]]))
    }
  }
  
  if(length(intG) == 0) {
    stop("The ID types of all omics data in the targetOmicsList 
         should be the same.")
  }
  
  #####process chrome location###
  
  chromeList <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                  "14","15","16","17","18","19","20","21","22","X","Y","All")
  
  x <- setdiff(chrome_sourceOmics, chromeList)
  if(length(x) > 0) {
    stop('The input chrome information for source omics data contains the 
         invalid information. Please only input the chromosome names from 
         the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
         "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
  }
  
  x <- setdiff(chrome_targetOmics, chromeList)
  if(length(x) > 0) {
    stop('The input chrome information for target omics data contains the 
         invalid information. Please only input the chromosome names from 
         the following list: "1","2","3","4","5","6","7","8","9","10","11","12",
         "13","14","15","16","17","18","19","20","21","22","X","Y" and "All".')
  }
  
  if((length(chrome_sourceOmics) > 1 && which(chrome_sourceOmics == "All") > 1) 
     || (length(chrome_sourceOmics) == 1 && chrome_sourceOmics == "All")) {
    chrome_sourceOmics <- "All"
  }
  
  if((length(chrome_targetOmics) > 1 && which(chrome_targetOmics == "All") > 1) 
     || (length(chrome_targetOmics) == 1 && chrome_targetOmics == "All")) {
    chrome_targetOmics <- "All"
  }
  
  if(chrome_sourceOmics == "All") {
    chrome_sourceOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                            "12","13","14","15","16","17","18","19","20","21","22","X","Y")
  }
  
  if(chrome_targetOmics == "All") {
    chrome_targetOmics <- c("1","2","3","4","5","6","7","8","9","10","11",
                            "12","13","14","15","16","17","18","19","20","21","22","X","Y")
  }
  
  #######Extract sub list#########
  genelocate <- datacache$genelocate
  
  genelocate_sourceOmics <- genelocate[genelocate[,2] %in% chrome_sourceOmics,]
  genelocate_targetOmics <- genelocate[genelocate[,2] %in% chrome_targetOmics,]
  
  intG <- intersect(intG, genelocate_targetOmics[,1])
  
  if(length(intG) == 0) {
    stop("The ID types for all omics data in the targetOmicsList should be 
         gene symbol or all genes in the target omics data are not in the 
         selected chromosomal location chrome_targetOmics.")
  }
  
  for(i in seq_len(length(targetOmicsList))) {
    targetOmicsList[[i]] <- targetOmicsList[[i]][intG,]
  }
  
  source_gene <- rownames(sourceOmics)
  source_gene_locate <- intersect(unique(genelocate_sourceOmics[,1]), source_gene)
  if(length(source_gene_locate) == 0) {
    stop("The ID type in the source omics data should be gene symbol or all 
         genes in the source omics data are not in the selected chromosomal 
         location chrome_sourceOmics.")
  }
  source_gene <- sourceOmics[source_gene_locate,]
  
  genelocate_sourceOmics <- genelocate_sourceOmics[genelocate_sourceOmics[,1] 
                                                   %in% source_gene_locate,]
  genelocate_targetOmics <- genelocate_targetOmics[genelocate_targetOmics[,1] 
                                                   %in% intG,]
  
  ###Calculate the correlation between cna and other omics data######
  cat("Identify the significant correlations...\n")
  if(length(targetOmicsList) == 1) {
    resultList <- calculateCorForTwoMatrices(source_gene, targetOmicsList[[1]], fdrThr)
  } else {
    cl <- makeCluster(nThreads)
    registerDoParallel(cl)
    resultList <- list()
    resultList <- foreach(i=seq_len(length(targetOmicsList)), 
                          .packages="multiOmicsViz") %dopar% {
                            corrArray <- calculateCorForTwoMatrices(source_gene, targetOmicsList[[i]], fdrThr)
                            return(corrArray)
                          }
    stopCluster(cl)
  }
  
  ##Calculate the location of genes in the heatmap
  chromLength <- datacache$chromLength
  
  re <- .calculateChromLength(chromLength, chrome_sourceOmics, genelocate_sourceOmics)
  genelocate_sourceOmics <- re$genelocate
  chromLength_sourceOmics <- re$chromLength
  
  re <- .calculateChromLength(chromLength, chrome_targetOmics, genelocate_targetOmics)
  genelocate_targetOmics <- re$genelocate
  chromLength_targetOmics <- re$chromLength
  
  ##########Plot Figure############
  
  cat("Plot figure...\n")
  
  if(length(targetOmicsList) == 1) {
    pdf(outputfile, height=8, width=8)
    .plotHeatMap(resultList, genelocate_sourceOmics, chromLength_sourceOmics,
                 genelocate_targetOmics, chromLength_targetOmics, sourceOmicsName,
                 targetOmicsName, dim=1)
    if(legend == TRUE) {
      legend("topleft", c("positive correlation", "negative correlation"),
             col=c("#FB6542", "#375E97
