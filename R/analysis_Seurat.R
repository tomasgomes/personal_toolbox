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
  obj = RunUMAP(obj, dims = 1:ncomp, verbose = F)
  obj = FindClusters(obj, algorithm = 2, verbose = F, graph.name = paste0(red, ncomp),
                     resolution = seq(0.2, 1.5, 0.1))
  # setting a more sensible identity as default
  obj = SetIdent(obj, value = paste0(red, ncomp, "_res.0.3"))
  return(obj)
}
