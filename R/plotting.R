##### FORMATTING #####
breakStr = function(s, n = 20) {return(gsub(paste0('(.{1,',as.character(n),'})(\\s|$)'), '\\1\n', s))}


##### QUALITY CONTROL #####

plotCountsGenes = function(obj){
  df = obj@meta.data[,c("nCount_RNA", "nFeature_RNA")]
  
  rr = diff(range(df$nCount_RNA))/diff(range(df$nFeature_RNA))
  cnt_g_plt = ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA))+
    scale_fill_gradient(trans = "log10")+
    labs(x = "# counts", y = "# features",
         fill = paste0("# cells\n(total: ", nrow(df), ")"))+
    theme_bw()+
    theme(aspect.ratio = 1,
          legend.title = element_text(size = 9),
          axis.text = element_text(colour = "black"))
  cnt_g_plt = if(nrow(df)<1500){
    cnt_g_plt+geom_point(size = 0.6, shape = 20, colour = "grey10", alpha = 0.2)
  } else{
    cnt_g_plt+geom_bin2d(binwidth = c(100*rr, 100))
  }
  
  rr = diff(range(log10(df$nCount_RNA)))/diff(range(log10(df$nFeature_RNA)))
  cnt_g_log = ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA))+
    scale_fill_gradient(trans = "log10")+
    scale_x_log10() + scale_y_log10()+
    labs(x = "# counts", y = "# features",
         fill = paste0("# cells\n(total: ", nrow(df), ")"))+
    theme_bw()+
    theme(aspect.ratio = 1,
          legend.title = element_text(size = 9),
          axis.text = element_text(colour = "black"))
  cnt_g_log = if(nrow(df)<25000){
    cnt_g_log+geom_point(size = 0.6, shape = 20, colour = "grey10", alpha = 0.2)
  } else{
    cnt_g_log+geom_bin2d(binwidth = c(100*rr, 100))
  }
  
  cnt_vln = ggplot(df, aes(x = "", y = nCount_RNA))+
    geom_violin(fill = "grey80")+ scale_y_log10()+
    geom_jitter(size = 0.3, alpha = 0.33)+
    labs(y = "# counts", x = "")+
    theme_bw()+
    theme(aspect.ratio = 3,
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank())
  
  fea_vln = ggplot(df, aes(x = "", y = nFeature_RNA))+
    geom_violin(fill = "grey80")+ scale_y_log10()+
    geom_jitter(size = 0.3, alpha = 0.33)+
    labs(y = "# features", x = "")+
    theme_bw()+
    theme(aspect.ratio = 3,
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank())
  
  res = list(cnt_g_plt = cnt_g_plt, cnt_g_log = cnt_g_log, cnt_vln = cnt_vln, fea_vln = fea_vln)
  return(res)
}

plotVar = function(obj, var.use = "percent.mt", lab.use = "MT genes [%]"){
  if(length(var.use)!=length(lab.use)){
    stop("var.use and lab.use must have the same length.")
  }
  res = list()
  for(v in 1:length(var.use)){
    vn = var.use[v]
    df = data.frame("v" = obj@meta.data[,vn])
    df = df[order(df$v, decreasing = F),,drop = F]
    
    plt = ggplot(df, aes(x = "", y = v))+
      geom_violin(fill = "grey80")+
      geom_jitter(size = 0.3, alpha = 0.33)+
      labs(y = lab.use[v], x = "")+
      theme_bw()+
      theme(aspect.ratio = 3,
            axis.text = element_text(colour = "black"),
            axis.ticks.x = element_blank())
    res[[vn]] = plt
  }

  return(res)
}



##### SPATIAL #####

# Improved FeaturePlot
# adapted from Beatriz Val
plot.exp <- function(obj, list.genes, sample_n = "", thr = 0){
  plot.list <- list()
  for (gene in list.genes){
    print(gene)
    plt <- SpatialFeaturePlot(obj, features = gene, alpha = 0.7, 
                              image.alpha = 0.8,   pt.size.factor = 2.5, 
                              keep.scale = 'all', images = "Slice") +
      scale_fill_viridis_c(option = 'A', end = 0.9, direction = -1, 
                           na.value = "#FFFFFF33") +
      theme(legend.position = "inside",
            legend.position.inside = c(1,0), legend.justification = c(1,0), 
            panel.grid = element_blank(), 
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8))
    
    if(sample_n!=""){
      plt=plt+labs(title = paste0(sample_n, ' ', gene))
    }
    
    plot.list[[gene]] <- plt
    
    # ensure that values below threshold (i.e. no expression) are NA to appear transparent
    plot.list[[gene]][[1]]$data[,gene][plot.list[[gene]][[1]]$data[,gene]<=thr] = NA
    plot.list[[gene]][[1]]$layers[[1]]$data[,gene][plot.list[[gene]][[1]]$layers[[1]]$data[,gene]<=thr] = NA
  }
  return(plot.list)
}


