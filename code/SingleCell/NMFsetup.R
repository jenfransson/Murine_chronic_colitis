
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


viridis_colors = c(viridis::viridis(7),viridis::viridis(8))
names(viridis_colors) = c(seq(0,60,10), 0:7)

showfactors_TCT = function(Groups, relvln = c(2.2,1)){
  factors = colnames(Groups[[1]]@meta.data)[
    grep("^nmf_", colnames(Groups[[1]]@meta.data))]
  
  Groups = lapply(Groups, function(x){
    AddMetaData(x, gsub("[WD]","",x$Time), "TimeNumber")
  })
  
  allmeta = do.call("rbind", lapply(Groups, function(x){x@meta.data}))
  
  groupnames = c("SCID","+ \U03B1  IL12","Rag2","+ \U03B1  IL12")
  names(groupnames) = c("SCID ut", "SCID il", "Rag2 ut", "Rag2 il")
  
  
  bars = lapply(as.numeric(gsub("nmf_","",factors)),function(i){
    barbreaks = seq(0,1,0.2)
    if(max(Groups[[1]]@reductions$NMF@feature.loadings[,i])<0.2){
      barbreaks = seq(0,1,0.1)
    }
    x = FactorGeneLoadingPlot(Groups[[1]], factor = i, topn = 10) +
      scale_y_continuous(expand = expansion(mult = 0),
                         breaks = barbreaks) + 
      theme(axis.text.y = element_text(face = "italic",
                                       size = 8, color = "black"),
            axis.text.x = element_text(size = 8, color = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8,hjust = 0),
            plot.margin = margin(t = 1, b = 1, l = 3, unit = "pt")) +
      labs(y = "Weight")
    x$layers[[1]]$aes_params$colour = NA
    x 
  })
  
  
  vlnlist = list()
  for(i in names(Groups)){
    message(i)
    y = VlnPlot(Groups[[i]],
                           features = factors,
                           group.by = "TimeNumber",
                           pt.size = 0.1,
                           combine = FALSE)
    y = lapply(y, function(x){
      x$layers[[2]]$aes_params$alpha =0.2
      x$layers[[2]]$aes_params$colour = "#888888"
      x$layers[[2]]$aes_params$shape = 16
      
      violin = x$layers[[1]]
      points = x$layers[[2]]
      violin$aes_params$alpha = 0.3
      violin$aes_params$size = 0.3
      x$layers[[1]] = rasterise(points, dpi = 300)
      x$layers[[2]] = violin
      
      if(max(x$data[,1])<1){
        ybreaks = c(0,0.4,0.8)
      }else{
        if(max(x$data[,1])<3.5){
          ybreaks = c(0:3)
        }else{
          if(max(x$data[,1])<7){
            ybreaks = seq(0,6,2)
          }else{
            if(max(x$data[,1])<10){
              ybreaks = seq(0,9,3)
            }else{
              if(max(x$data[,1])<16){
                ybreaks = seq(0,15,5)
              }else{
                ybreaks = seq(0,100,10)
              }
            }
          }
        }
      }
      
      x + NoLegend() + scale_y_continuous(breaks = ybreaks) + 
        theme(axis.title.y = element_blank(),
              axis.title.x = element_text(
                size = 8,
                margin = margin(t = 1, r = 1, b = 1, l = 1,
                                unit = "pt"),
                hjust = 0),
              plot.title = element_text(
                size = 8, hjust = 0,
                margin = margin(t = 3, r = 1, b = 1, l = 1, unit = "pt"),
                face = "bold"
              ),
              axis.text = element_text(size = 8),
              axis.text.x = element_text(angle = 0, 
                                         hjust = 0.5)) +
        scale_fill_manual(values = viridis_colors)
    })
    
    vlnlist[[i]] = y
  }
  
  gridlist = list()
  
  for(i in 1:length(vlnlist[[1]])){
    message(i)
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
        
        colnames(df)[1] <- "y"
        df.summary <- df %>%
          group_by(ident) %>%
          summarise(
            sd = sd(y, na.rm = TRUE),
            m = mean(y)
          )
        
        
        
        x = x + 
          
          labs(x = paste(c("Weeks","Days")[model],
                         "post naive CD4+ TCT")) +
          coord_cartesian(ylim=c(0,ymax*1.1)) +
          theme(plot.margin = margin(t = c(0,2)[model], r = 1, b = c(2,0)[model], l = 1, unit = "pt")) + 
          geom_point(aes(y = m),data = df.summary, 
                     size = 1) + 
          labs(title = groupnames[[n]])
        if(n %in% c("Rag2 il","SCID il")){
          x = x + theme(axis.text.y = element_blank(),
                        axis.title.x = element_blank())
        }
        x
      }), rel_widths = relvln, align = "h"), ncol = 1#,
    ))
  }
  
  plots <- lapply(1:length(gridlist), function(i){
    print(plot_grid(
      ggdraw() + 
        draw_label(
          "Factor score",
          hjust = 0.5, 
          size =  8, angle = 90
        ),
      plot_grid(ggdraw() + 
                  draw_label(
                    i,
                    fontface = 'bold',
                    hjust = 0.5, 
                    size =  10
                  ),
                plot_grid(gridlist[[i]] + theme(plot.title = element_blank()),
                          bars[[i]], ncol = 2, rel_widths = c(2,1)), 
                rel_heights = c(1,6), ncol = 1), nrow = 1, rel_widths = c(1,10)))
  })
}


showfactors_Il10 = function(nmf){
  factors = colnames(nmf@meta.data)[grep("^nmf_", colnames(nmf@meta.data))]
  vlns = lapply(VlnPlot(nmf,
                        features = factors,
                        group.by = "G_W",
                        pt.size = 0.1,
                        combine = FALSE),
                function(x){
                  df = x$data
                  df$Week = factor(as.numeric(gsub('.* ','',df$ident)))
                  x$data = df
                  x$layers[[1]]$mapping$fill <- as.name("Week")
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
                  
                  x = x + NoLegend() +
                    labs(x = "Week") +
                    theme(axis.title.y = element_blank(),
                          axis.title.x = element_text(size = 8,
                                                      margin = margin(t = 1, r = 1, b = 1, l = 1,
                                                                      unit = "pt"),
                                                      hjust = 0.5),
                          plot.title = element_text(size = 8, hjust = 0.5,
                                                    margin = margin(t = 3, r = 1, b = 1, l = 1, unit = "pt"),
                                                    face = "bold"
                          ),
                          #plot.title = element_blank(),
                          axis.text = element_text(size = 8),
                          axis.text.x = element_text(angle = 0, 
                                                     hjust = 0.5),
                          plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "pt")) + 
                    #geom_errorbar(aes(ymin = m, ymax = m+sd, y = m),data = df.summary, 
                    #          #      color = "red", 
                    #              width = 0.25) +
                    geom_point(aes(y = m),data = df.summary, 
                               #     color = "red", 
                               size = 1) + 
                    scale_x_discrete(labels = c("6 (WT)",6,8,10,12,14)) +
                    scale_fill_viridis_d()
                  
                  x
                })
  
  bars = lapply(as.numeric(gsub("nmf_","",factors)),function(i){
    barbreaks = seq(0,1,0.2)
    if(max(nmf@reductions$NMF@feature.loadings[,i])<0.2){
      barbreaks = seq(0,1,0.1)
    }
    x = FactorGeneLoadingPlot(nmf, factor = i, topn = 10) +
      scale_y_continuous(expand = expansion(mult = 0),
                         breaks = barbreaks) + 
      theme(axis.text.y = element_text(face = "italic",
                                       size = 8, color = "black"),
            axis.text.x = element_text(size = 8, color = "black"),
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 8,hjust = 0),
            plot.margin = margin(t = 1, b = 1, l = 3, unit = "pt")) +
      labs(y = "Weight")
    x$layers[[1]]$aes_params$colour = NA
    #x$layers[[1]]$aes_params$fill = "#888888"
    x 
  })
  
  plots <- lapply(1:length(vlns), function(i){
    #lapply(1, function(i){
    print(plot_grid(
      ggdraw() + 
        draw_label(
          "Factor score",
          #fontface = 'bold',
          hjust = 0.5, 
          size =  8, angle = 90
        ),
      plot_grid(ggdraw() + 
                  draw_label(
                    i,
                    fontface = 'bold',
                    hjust = 0.5, 
                    size =  10
                  ),
                plot_grid(vlns[[i]] + theme(plot.title = element_blank()),
                          bars[[i]], ncol = 2, rel_widths = c(1.5,1)), 
                rel_heights = c(1,6), ncol = 1), nrow = 1, rel_widths = c(1,10)))
  })
  
  
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
