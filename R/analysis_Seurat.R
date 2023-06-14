##### STANDARD ANALYSIS #####

basicSeurat = function(obj){
  obj = suppressWarnings(SCTransform(obj, do.correct.umi = T, vars.to.regress = "nCount_RNA",
                                     verbose = F, variable.features.rv.th = 1,
                                     variable.features.n = NULL, seed.use = 1))
  obj = RunPCA(obj, verbose = F)
  return(obj)
}

runSeuratClust = function(obj, red = "pca", ncomp = 10){
  if(ncomp=="all"){ncomp = ncol(Embeddings(obj, red))}
  obj = FindNeighbors(obj, dims = 1:ncomp, force.recalc = T, verbose = F,
                      reduction = red, graph.name = paste0(red, ncomp))
  obj = RunUMAP(obj, dims = 1:ncomp, reduction = red, verbose = F)
  obj = FindClusters(obj, algorithm = 2, verbose = F, graph.name = paste0(red, ncomp),
                     resolution = seq(0.2, 1.5, 0.1))
  # setting a more sensible identity as default
  obj = SetIdent(obj, value = paste0(red, ncomp, "_res.0.3"))
  return(obj)
}


# remove normalisation (go from data to counts in RNA assay)
## works cell by cell
removeNorm = function(r, N = 10) {
  r = expm1(r)
  h = as.data.frame(table(r))
  h$r = as.numeric(as.character(h$r))
  h = h[order(h$r, decreasing = F),]
  N_ = min(nrow(h), N)
  y = head(h$r, N_)
  x = seq_len(N_)
  s = coef(lm(x ~ y))[2]
  
  return(r*s)
}



# Replace gene names in different slots of a Seurat object.
## gene_match_table has two columns: original_names and new_names
RenameGenesSeurat = function(obj, gene_match_table){ 
  # change assays
  for(n in SeuratObject::Assays(obj)){
    g_use = gene_match_table[match(rownames(obj[[n]]@counts), gene_match_table[,1]),2]
    rownames(obj[[n]]@counts) = g_use
    g_use = gene_match_table[match(rownames(obj[[n]]@data), gene_match_table[,1]),2]
    rownames(obj[[n]]@data) = g_use
    g_use = gene_match_table[match(rownames(obj[[n]]@scale.data), gene_match_table[,1]),2]
    rownames(obj[[n]]@scale.data) = g_use
    g_use = gene_match_table[match(rownames(obj[[n]]@meta.features), gene_match_table[,1]),2]
    rownames(obj[[n]]@meta.features) = g_use
    g_use = gene_match_table[match(obj[[n]]@var.features, gene_match_table[,1]),2]
    obj[[n]]@var.features = g_use
  }
  
  # change reductions
  for(n in SeuratObject::Reductions(obj)){
    g_use = gene_match_table[match(rownames(obj[[n]]@feature.loadings), gene_match_table[,1]),2]
    rownames(obj[[n]]@feature.loadings) = g_use
    g_use = gene_match_table[match(rownames(obj[[n]]@feature.loadings.projected), 
                                   gene_match_table[,1]),2]
    rownames(obj[[n]]@feature.loadings.projected) = g_use
  }
  
  return(obj)
}



# Make prompt to ChatGPT
# presto_df is the data frame resulting from running presto
# add_info is a named list to input additional info (see example)
# example annotationPrompt(deg, n = 10, add_info = list("species" = "human", "tissue" = "liver"))
annotationPrompt = function(presto_df, n = 10, 
                            pval_thr = 0.05, pct_diff = 20, logFC_thr = 0.2,
                            add_info = NULL){
  require(dplyr)
  
  topmk = presto_df |>
    group_by(group) |> # for each cluster
    filter(padj<=pval_thr) |> # only p-value below 0.05
    filter(pct_in>(pct_out+pct_diff)) |> 
    filter(logFC>logFC_thr) |> 
    top_n(n = n, wt = logFC) # top genes per cluster, ranked by logFC
  
  prompt = "I have performed a scRNA-seq analysis, and have encountered various clusters, for which I calculate the marker genes. The top marker genes for each cluster are as follows:\n"
  for(cl in unique(topmk$group)){
    prompt = paste0(prompt, " - cluster ", cl, ": ", 
                    paste0(topmk$feature[topmk$group==cl], collapse = ", "), "\n")
  }
  
  prompt = paste0(prompt, 
                  "\nHere is more information about the data: ", 
                  paste0(sapply(names(add_info), function(x) paste0(x, ": ", add_info[x])), 
                         collapse = "; "), "\n")
  
  prompt = paste0(prompt, 
                  "Can you tell me what are the most likely cell types that each cluster matches to?")
  cat(prompt)
}



