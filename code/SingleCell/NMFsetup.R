
suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
  library(dplyr)
  library(Matrix)
  library(gprofiler2)
  library(inflection)
  library(viridis)
  library(ggrastr)
})

source("RunNMF.R")
source("gaussian_smoothing.R")

set.seed(1)


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
viridis_colors = unique(c(viridis::viridis(7),viridis::viridis(8)))
names(viridis_colors) = unique(c(as.character(seq(0,1,length.out = 7)),
                                 as.character(seq(0,1,length.out = 8))))

makeViolins_Tcelltransfer <- function(
  Groups, features, relvln = c(2.2,1), ybreaks = 0:10){
  
  allmeta = do.call("rbind", lapply(Groups, function(x){x@meta.data}))
  
  groupnames = c("SCID","+ \U03B1  IL12","Rag2","+ \U03B1  IL12")
  names(groupnames) = c("SCID ut", "SCID il", "Rag2 ut", "Rag2 il")
  
  vlnlist = list()
  for(i in names(Groups)){
    message(i)
    vlnlist[[i]] = VlnPlot(Groups[[i]],
                           features = features,
                           group.by = "TimeNumber",
                           pt.size = 0.1,
                           combine = FALSE)
  }
  
  gridlist = list()
  
  for(i in 1:length(vlnlist[[1]])){
    plots =  lapply(vlnlist, `[[`, i)
    names(plots) = names(vlnlist)
    ymax = max(sapply(plots, function(x){
      max(x$data[,1])
    }))
    
    title <- ggdraw() + 
      draw_label(
        plots[[1]]$labels$title,
        fontface = 'bold',
        hjust = 0.5, 
        size =  10
      )
    
    gridlist[[i]] = (plot_grid(#title,
      plot_grid(plotlist = lapply(names(plots), function(n){
        x = plots[[n]]
        model = (substr(n,1,4) == "Rag2")+1
        df = x$data
        df$RelTime <- as.character(allmeta[rownames(df),"RelTime"])
        x$data = df
        x$layers[[1]]$mapping$fill <- as.name("RelTime")
        x$mapping$fill = NULL
        
        colnames(df)[1] <- "y"
        df.summary <- df %>%
          group_by(ident) %>%
          summarise(
            sd = sd(y, na.rm = TRUE),
            m = mean(y)
          )
        
        x$layers[[2]]$aes_params$alpha =0.2
        x$layers[[2]]$aes_params$colour = "#888888"
        x$layers[[2]]$aes_params$shape = 16
        
        
        
        violin = x$layers[[1]]
        points = x$layers[[2]]
        #points$aes_params$alpha = 0.4
        violin$aes_params$alpha = 0.3
        violin$aes_params$size = 0.3
        x$layers[[1]] = rasterise(points, dpi = 300)
        x$layers[[2]] = violin
        
        x = x + NoLegend() + scale_y_continuous(breaks = ybreaks) + 
          coord_cartesian(ylim=c(0,ymax*1.1)) +
          labs(x = paste(c("Weeks","Days")[model],
                         "post naive CD4+ TCT")) +
          theme(axis.title.y = element_blank(),
                axis.title.x = element_text(size = 8,
                                            margin = margin(t = 1, r = 1, b = 1, l = 1,
                                                            unit = "pt"),
                                            hjust = 0),
                plot.title = element_text(size = 8, hjust = 0,
                                          margin = margin(t = 3, r = 1, b = 1, l = 1, unit = "pt"),
                                          face = "bold"
                ),
                #plot.title = element_blank(),
                axis.text = element_text(size = 8),
                axis.text.x = element_text(angle = 0, 
                                           hjust = 0.5),
                plot.margin = margin(t = c(0,2)[model], r = 1, b = c(2,0)[model], l = 1, unit = "pt")) + 
          #geom_errorbar(aes(ymin = m, ymax = m+sd, y = m),data = df.summary, 
          #          #      color = "red", 
          #              width = 0.25) +
          geom_point(aes(y = m),data = df.summary, 
                     #     color = "red", 
                     size = 1) + 
          labs(title = groupnames[[n]]) +
          scale_fill_manual(values = viridis_colors)
        if(n %in% c("Rag2 il","SCID il")){
          x = x + theme(axis.text.y = element_blank(),
                        axis.title.x = element_blank())
        }
        x
      }), rel_widths = relvln, align = "h"), ncol = 1#,
      #rel_heights = c(1,10)
    ))
  }
  
  return(gridlist)
}



makeplots = function(pwys){
  lapply(1:ntop, function(i) {
    if (!paste0("factor_", i) %in% pwys$factor) return(NULL)
    topdata = subset(pwys, 
                     factor %in% paste0("factor_", i))
    if(diff(range(topdata$GeneRatio))<0.1){
      sizebreaks = waiver()
    }else{
      if(diff(range(topdata$GeneRatio))>0.4){
        sizebreaks = seq(0,1,0.2)
      }else{
        sizebreaks = seq(0,1,0.1)
      }}
    txtoffset = diff(range(-log10(topdata$p_value)))/10
    if(nrow(topdata)>10){topdata = topdata[1:10,]}
    g <- ggplot(data = topdata, 
                aes(reorder(term_name, -log10(p_value)), -log10(p_value),
                    fill = -log10(p_value))) +
      geom_point(aes(size = GeneRatio),
                 color = "black", shape = 21) +
      geom_text(mapping = aes(label = paste0(intersection_size),#,"/",query_size), 
                              y = -log10(p_value) + txtoffset),
                hjust = 0,
                show.legend = FALSE, size = 3) + 
      coord_flip(clip = "off") +
      #facet_grid(~factor) +
      scale_size_continuous(range = c(0.5, 6),
                            breaks = sizebreaks,
                            limits = c(0,NA)) +
      scale_fill_gradientn(colours = viridis::magma(n = 9) %>% rev(),
                           breaks = scales::trans_breaks(identity, identity, n= 3)) +
      theme_minimal() +
      theme(legend.position = "bottom", 
            axis.title.y = element_blank(),
            axis.text = element_text(size = 8),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      guides(
        size = guide_legend(title.position = "top", title.hjust = 0.5,
                            title.theme = element_text(size = 8),
                            text.theme = element_text(size = 8)),
        fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                              title.theme = element_text(size = 8),
                              text.theme = element_text(size = 8),
                              barheight = unit(6,"pt"),
                              barwidth = unit(50,"pt"))) +
      labs(x = "term", y = "", title = paste0("ent_f",i," (",topdata$query_size[[1]]," genes)"))
    g
  })
}

shortenterms = function(terms){
  sapply(strsplit(terms," "), function(term){
    concat = term[[1]]
    if(length(term)>1){
      currentlength = nchar(concat)
      for(i in 2:length(term)){
        if(currentlength>25){
          concat = paste(concat, term[[i]], sep="\n")
          currentlength = 0
        }else{
          concat = paste(concat, term[[i]])
          currentlength = currentlength + nchar(term[[i]]) + 1
        }
        
        
      }
    }
    concat
  })
}



getclusternames_Il10 = function(celltype, clusters){
  stop("Get cluster names in file (once)")
  
  # clusternames = read.csv2("../Clusternames_IL10KO.csv")
  # #clusternames = data.frame(RNA_clusters = 0:20, epi_cluster = 0:20, epi_order = 0:20)
  # 
  # clusterorder_all = clusternames[, paste0(celltype,"_order")]
  # clusterorder = clusterorder_all[clusterorder_all %in% 
  #                                   clusternames[clusternames$Cluster %in% clusters, 
  #                                 paste0(celltype,"_clusters")]]
  # 
  # factor(clusternames[match(clusters, clusternames[,1]),
  #                     paste0(celltype,"_clusters")], 
  #        levels = clusterorder)
  
}

getfactorpvalues_TCT = function(nmf, downsample_n, ntop){ ### Added 2024-01-02
  set.seed(42)
  nmf = SetIdent(nmf, value = "orig.ident")
  subsample = subset(nmf, downsample = downsample_n)
  
  subsample.metadata = subsample@meta.data
  
  topicnumbers = c(1:ntop)
  names(topicnumbers) = paste0("nmf_", topicnumbers)
  modelpvalues = as.data.frame(lapply(topicnumbers, function(x){
    y = (summary(lm(as.formula(paste0(paste0("nmf_",x,"~"), "Model*RelTime")),
                    data = subsample.metadata[subsample$Treatment!="il",]))$coefficients[-1,4]);
    z = as.numeric(y); names(z) = names(y);z}))
  modelpvalues_corr = adjustvalues(modelpvalues)
  modelpvalues_corr$depvar = rownames(modelpvalues_corr)
  
  return(modelpvalues_corr)
}


getfactorpvalues_Il10 = function(nmf, downsample_n, ntop){ ### Added 2024-01-02
  set.seed(1234)
  nmf = SetIdent(nmf, value = "orig.ident")
  subsample = subset(nmf, downsample = downsample_n)
  
  subsample.metadata = subsample@meta.data
  
  topicnumbers = c(1:ntop)
  names(topicnumbers) = paste0("nmf_", topicnumbers)
  modelpvalues = as.data.frame(lapply(topicnumbers, function(x){
    y = (summary(lm(as.formula(paste0(paste0("nmf_",x,"~"), "Week")),
                    data = subsample.metadata[
                      subsample$Genotype!="WT",]))$coefficients[
                        -1,4, drop = FALSE]);
    z = as.numeric(y); names(z) = rownames(y);z}))
  modelpvalues_corr = adjustvalues(modelpvalues)
  modelpvalues_corr$depvar = rownames(modelpvalues_corr)
  
  rank_modelpvalues = as.data.frame(lapply(topicnumbers, function(x){
    y = (summary(lm(as.formula(paste0(paste0("rank(nmf_",x,")~"), "Week")),
                    data = subsample.metadata[
                      subsample$Genotype!="WT",]))$coefficients[
                        -1,4, drop = FALSE]);
    z = as.numeric(y); names(z) = rownames(y);z}))
  
  rank_modelpvalues_corr = adjustvalues(rank_modelpvalues)
  rank_modelpvalues$depvar_corr = rownames(rank_modelpvalues_corr)
  
  return(list(modelpvalues = modelpvalues_corr, 
              rank_modelpvalues = rank_modelpvalues))
}

adjustvalues = function(testframe){
  as.data.frame(matrix(p.adjust(as.matrix(testframe),"BH"), 
                       nrow = nrow(testframe), 
                       dimnames = list(rownames(testframe),
                                       colnames(testframe))))
}