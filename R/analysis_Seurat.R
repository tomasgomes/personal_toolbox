##### STANDARD ANALYSIS #####

# obtain basic Seurat QC metrics
## setub mostly for human, and genes in ENS000-NAME format
basicQCSeurat = function(obj, mt.pat = "-MT-", calc.ribo = TRUE, g.pat = c("","")){
  # ribosomal genes
  ## obtained from https://www.genenames.org/data/genegroup/#!/group/728
  ## and https://www.genenames.org/data/genegroup/#!/group/729
  rb_genes = c("RPLP0","RPLP1","RPLP2","RPL3","RPL3L","RPL4","RPL5","RPL6",
               "RPL7","RPL7A","RPL7L1","RPL8","RPL9","RPL10","RPL10A","RPL10L","RPL11",
               "RPL12","RPL13","RPL13A","RPL14","RPL15","RPL17","RPL18","RPL18A","RPL19",
               "RPL21","RPL22","RPL22L1","RPL23","RPL23A","RPL24","RPL26","RPL26L1","RPL27",
               "RPL27A","RPL28","RPL29","RPL30","RPL31","RPL32","RPL34","RPL35","RPL35A",
               "RPL36","RPL36A","RPL36AL","RPL37","RPL37A","RPL38","RPL39","RPL39L","UBA52",
               "RPL41","RPSA","RPS2","RPS3","RPS3A","RPS4X","RPS4Y1","RPS4Y2","RPS5","RPS6",
               "RPS7","RPS8","RPS9","RPS10","RPS11","RPS12","RPS13","RPS14","RPS15",
               "RPS15A","RPS16","RPS17","RPS18","RPS19","RPS20","RPS21","RPS23","RPS24",
               "RPS25","RPS26","RPS27","RPS27A","RPS27L","RPS28","RPS29","FAU")
  rb_genes_obj = unlist(sapply(rb_genes, 
                               function(x) grep(paste0(g.pat[1],x,g.pat[2]), 
                                                rownames(obj[["RNA"]]), 
                                                value = T)))
  
  # scores
  if(!isnull(mt.pat)){
    obj = Seurat::PercentageFeatureSet(obj, pattern = mt.pat, 
                                       col.name = "percent.mt", assay = "RNA")
  }
  if(calc.ribo){
    obj = Seurat::PercentageFeatureSet(obj, col.name = "percent.rb", 
                                       assay = "RNA", features = rb_genes_obj)
  }
  
  # cell cycle scoring
  # might have trouble with bins so running iteratively until it works
  s.genes = paste0(g.pat[1], cc.genes$s.genes, g.pat[2])
  s.genes_obj = unlist(sapply(s.genes, 
                              function(x) grep(x, rownames(obj[["RNA"]]), 
                                               value = T)))
  g2m.genes = paste0(g.pat[1], cc.genes$g2m.genes, g.pat[2])
  g2m.genes_obj = unlist(sapply(g2m.genes, 
                                function(x) grep(x, rownames(obj[["RNA"]]), 
                                                 value = T)))
  bins = 24
  while(T){
    objwsc = tryCatch({Seurat::CellCycleScoring(obj, nbin = bins, assay = "RNA",
                                                s.features = s.genes_obj,
                                                g2m.features = g2m.genes_obj,
                                                set.ident = F)},
                      error = function(e){print(bins)})
    
    if(is.numeric(objwsc) & bins>6){
      bins = bins-2
    } else{
      obj = objwsc
      break
    }
  }
  
  return(obj)
}


# apply QC filters
applyQC = function(obj, thr_list){
  cond = rep(T, ncol(obj))
  l = list()
  for(n in names(thr_list)){
    l[[n]] = (obj@meta.data[,n]>thr_list[[n]][1]) & 
      (obj@meta.data[,n]<thr_list[[n]][2])
    cond = cond & 
      (obj@meta.data[,n]>thr_list[[n]][1]) & 
      (obj@meta.data[,n]<thr_list[[n]][2])
  }
  print(sapply(l, table))
  
  message(paste0("Removed ", sum(!cond), " out of ", length(cond), " cells."))
  
  return(obj[,cond])
}


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


# use DoubletFinder to classify droplets as doublets
removeDoublets = function(obj, ncomp = 10, sct = T, filter.obj = T){
  require(DoubletFinder)
  
  # parameter sweep for optimal pK
  sweep.res = paramSweep_v3(obj, PCs = 1:ncomp, sct = sct)
  sweep.stats = summarizeSweep(sweep.res, GT = FALSE)
  bcmvn = find.pK(sweep.stats)
  bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk = bcmvn.max$pK
  optimal.pk = as.numeric(levels(optimal.pk))[optimal.pk]
  
  # homeotypic doublets
  homotypic.prop = modelHomotypic(obj$seurat_clusters)
  nExp_poi = round(homotypic.prop*ncol(obj))
  nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
  
  #classify
  res = doubletFinder_v3(obj, PCs = 1:ncomp, pK = optimal.pk, 
                         nExp = nExp_poi, 
                         reuse.pANN = FALSE, sct = sct)
  
  if(filter.obj){
    varuse = grep("DF.classifications_", colnames(res@meta.data), value = T)
    res = res[,res@meta.data[,varuse]=="Singlet"]
  }
  
  return(res)
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


