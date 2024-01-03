# Make prompt to ChatGPT
# presto_df is the data frame resulting from running presto
# add_info is a named list to input additional info (see example)
# example annotationPrompt(deg, n = 10, add_info = list("species" = "human", "tissue" = "liver"))
annotationPrompt = function(presto_df, n = 10, 
                            pval_thr = 0.05, pct_diff = 20, logFC_thr = 0.2,
                            add_info = NULL, print.it = T){
  require(dplyr)
  
  topmk = presto_df |>
    group_by(group) |> # for each cluster
    filter(padj<=pval_thr) |> # only p-value below 0.05
    filter(pct_in>(pct_out+pct_diff)) |> 
    filter(logFC>logFC_thr) |> 
    arrange(logFC) |>
    top_n(n = n, wt = logFC) # top genes per cluster, ranked by logFC
  
  prompt = "I have performed a scRNA-seq analysis, and have encountered various clusters, for which I calculated the marker genes. The top marker genes for each cluster are as follows:\n"
  for(cl in unique(topmk$group)){
    prompt = paste0(prompt, " - cluster ", cl, ": ", 
                    paste0(topmk$feature[topmk$group==cl], collapse = ", "), "\n")
  }
  
  if(!is.null(add_info)){
    prompt = paste0(prompt, 
                    "\nTo interpret these clusters, these additional sample characteristics should be considered: ", 
                    paste0(sapply(names(add_info), function(x) paste0(x, ": ", add_info[x])), 
                           collapse = "; "), "\n")
  }
  
  prompt = paste0(prompt, 
                  "Can you tell me what are the most likely cell types that each cluster matches to? And report the top 5 genes used to reach this annotation.")
  
  if(print.it){
    cat(prompt)
  }
  return(prompt)
}


# calculate enriched genes in another set (e.g. GO Term in markers)
## written by ChatGPT
hypergeometricTest <- function(geneSet, backgroundSet, universeSize) {
  # Calculate the number of genes in the intersection of geneSet and backgroundSet
  intersectionSize <- length(intersect(geneSet, backgroundSet))
  
  # Calculate the number of genes in geneSet
  geneSetSize <- length(geneSet)
  
  # Calculate the number of genes in backgroundSet
  backgroundSetSize <- length(backgroundSet)
  
  # Calculate the number of genes not in either geneSet or backgroundSet
  notInSetSize <- universeSize - (geneSetSize + backgroundSetSize - intersectionSize)
  
  # Perform the hypergeometric test
  p_value <- phyper(intersectionSize - 1, geneSetSize, notInSetSize, backgroundSetSize, 
                    lower.tail = FALSE)
  
  # Create a result object with relevant information
  result <- list(
    geneSetSize = geneSetSize,
    backgroundSetSize = backgroundSetSize,
    intersectionSize = intersectionSize,
    notInSetSize = notInSetSize,
    p_value = p_value
  )
  
  return(result)
}


