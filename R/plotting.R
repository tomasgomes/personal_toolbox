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
    geom_jitter(size = 0.33, alpha = 0.33)+
    labs(y = "# counts", x = "")+
    theme_bw()+
    theme(aspect.ratio = 3,
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank())
  
  fea_vln = ggplot(df, aes(x = "", y = nFeature_RNA))+
    geom_violin(fill = "grey80")+ scale_y_log10()+
    geom_jitter(size = 0.33, alpha = 0.33)+
    labs(y = "# features", x = "")+
    theme_bw()+
    theme(aspect.ratio = 3,
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank())
  
  res = list(cnt_g_plt = cnt_g_plt, cnt_g_log = cnt_g_log, cnt_vln = cnt_vln, fea_vln = fea_vln)
  return(res)
}

plotVar = function(obj, var.use = "percent.mt", lab.use = "MT genes [%]"){
  df = data.frame("percent.mt" = obj@meta.data[,var.use])
  df = df[order(df$percent.mt, decreasing = F),,drop = F]
  
  mt_plt = ggplot(df, aes(x = "", y = percent.mt))+
    geom_violin(fill = "grey80")+
    geom_jitter(size = 0.33, alpha = 0.33)+
    labs(y = lab.use, x = "")+
    theme_bw()+
    theme(aspect.ratio = 3,
          axis.text = element_text(colour = "black"),
          axis.ticks.x = element_blank())
  
  res = list(mt_plt = mt_plt)
  return(res)
}


