
#### Shared for bulk ####


getgtf = function(datadir){
  if(!file.exists( paste0(datadir,"../Mus_musculus.GRCm38.101.gtf") )){
    scriptwd = getwd()
    setwd(paste0(datadir,"/.."))
    system("wget ftp://ftp.ensembl.org/pub/release-101/gtf/mus_musculus/Mus_musculus.GRCm38.101.gtf.gz")
    system("gunzip Mus_musculus.GRCm38.101.gtf.gz")
    setwd(scriptwd)
  }
  
  gtf <- import.gff(paste0(datadir,"../Mus_musculus.GRCm38.101.gtf") )
  gtf[ ! is.na(gtf$transcript_id) ,]
}

shortenterms = function(term){
  abbr = list(c("positive","pos"),
       c("negative","neg"),
       c("regulation","reg"),
       c("dependent","dep"))
  for(ab in abbr){
    term = gsub(ab[1],ab[2],term)
  }
  term
}

# Visualization parameters

{
treatmentshapes = c(21,22)
plottheme = theme_classic() + 
  theme(axis.title = element_text(size = 8),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.margin = margin(5.5,5.5,5.5,5.5,"pt"),
        legend.margin = margin(0, 0, 0, 0, "cm"), 
        legend.box.margin = margin(0, 0, 0, 0, "cm"), 
        axis.line = element_line(linewidth = 0.25, color = "black"))
pointsize = 1.8
strokewidth = 0.3
viridistype = "viridis"
viridisdirection = 1
treatmentcolors = c("#FF0000", "#000000")

bartheme = theme(panel.background = element_blank(), panel.grid = element_blank(),
                 axis.text.y = element_blank(), axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(), legend.position = "none",
                 axis.text.x = element_text(size = 6),axis.title.x = element_text(size = 6),
                 plot.title = element_text(size = 8)
)




linetheme_human = theme(axis.title = element_text(size = 8),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(size = 8),
                  axis.ticks = element_blank(),
                  panel.background = element_blank(),
                  plot.title = element_text(size = 10))

linetheme_aIL12 = theme(axis.title.x = element_text(size = 8),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.text.x = element_text(size = 8),
                        axis.ticks = element_blank(),
                        panel.background = element_blank(),
                        legend.position = "none",
                        legend.key = element_blank(),
                        legend.title = element_text(size = 8),
                        legend.text = element_text(size = 8))

bartheme_human = theme(axis.title = element_text(size = 8),
                 axis.text = element_text(size = 8),
                 axis.ticks.x = element_blank(),
                 plot.title = element_text(size = 10, hjust = 0.5),
                 panel.background = element_blank()
)
}


# gg_color_hue returns a list of n colors that are taken from the same 
#    palette as used by ggplot

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


# loadGeneSets creates or loads .rds-files with KEGG pathways and GO terms.
#    It returns a list containing two lists named "KEGG" and "BP". Each list contains
#    vectors of genes including in a term/pathway and with a name corresponding to the
#    term/pathway. These lists are formatted to be used as "pathways" in ownORA.

loadKEGGGeneSets = function(filepath_kegg = "mouseKEGG.rds",
                            savePathways = TRUE){
  if(file.exists(filepath_kegg) & savePathways){
    kegg_pathways = readRDS(filepath_kegg)
  }else{
    
    library(KEGGREST)
    
    mmupaths <- keggList("pathway",organism = "mmu")
    
    
    mmupathsinfo = lapply(names(mmupaths), function(x){
      res = keggGet(x)[[1]]
      return(list(res$GENE,res$NAME))
    })
    
    mmupathgenesclean = lapply(mmupathsinfo, function(i){
      if(!is.null(i[[1]])){
        geneinput <- i[[1]][seq(2,length(i[[1]]),2)]
        return(gsub(";.*", "", geneinput))
      }else{
        return(NULL)
      }
    })
    names(mmupathgenesclean) <- lapply(mmupathsinfo,"[[",2)
    
    kegg_pathways <- mmupathgenesclean
    
    if(savePathways) saveRDS(kegg_pathways, file = filepath_kegg)
  }
  
  return(kegg_pathways)
}

loadGOGeneSets = function(filepath_gobp = "mouseGOBP.rds",
                            savePathways = TRUE,
                            allgenes,
                            min_go_thresh = 3,
                            max_go_thresh = 3000){

  if(file.exists(filepath_gobp) & savePathways){
    gobp_pathways = readRDS(file = filepath_gobp)
  }else{
    
    library(org.Mm.eg.db)
    
    gobp_pathways_entrez <- list()
    
    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    goterms <- goterms[goterms == "BP"]
    
    go2genes <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, keys=names(goterms),
                                                       column="SYMBOL",
                                                       keytype="GOALL", multiVals='list'))
    
    gobp_pathways_all <- go2genes
    names(gobp_pathways_all) <- Term(GO.db::GOTERM)[match(names(go2genes), keys(GO.db::GOTERM))]
    
    presentgenes <- sapply(gobp_pathways_all, function(x){sum(!is.na(match(allgenes, x)))})
    gobp_pathways <- gobp_pathways_all[presentgenes>min_go_thresh & 
                                         presentgenes<max_go_thresh]
    
    if(savePathways) saveRDS(gobp_pathways, file = filepath_gobp)
    
  }
  
  return(gobp_pathways)
}



# ownORA performs an over-representation analysis based on a vector of genes, 
#    a list of pathways and a vector of background genes. It returns a data.frame
#    with the results for each pathway.

# genelist: A vector of gene names of interest. Names should match gene names in pathways list entries.
# pathways: A list containing groups of genes. Each list item should be a vector of genes, with the 
#    name of each entry being the name of the pathway. Lists of KEGG pathways and GO biological process
#    can be created using the loadGeneSets function.
# universe: A vector of gene names for background genes to test against. This could be e.g. all expressed 
#    genes in the dataset.

ownORA <- function(genelist, pathways, universe){
  
  ### Adapted from enricher_internal from DOSE
  
  intersects <- list()
  for(i in 1:length(pathways)){
    intersects[[i]] <- list(intersect(pathways[[i]], genelist),intersect(pathways[[i]], universe))
  }
  names(intersects) <- names(pathways)
  args.df <- data.frame(numWdrawn = sapply(intersects,function(x){length(x[[1]])})-1, ## White balls drawn
                        numW = sapply(intersects,function(x){length(x[[2]])}),        ## White balls
                        numB = length(unique(universe))-sapply(intersects,function(x){length(x[[2]])}),               ## Black balls
                        numDrawn = length(genelist))                        ## balls drawn
  rownames(args.df) <- names(pathways)
  args.df <- args.df[args.df[,2]>3,]
  
  pvalues <- apply(args.df, 1, function(n)
    phyper(n[1], n[2], n[3], n[4], lower.tail=FALSE)
  )
  terms = rownames(args.df)
  overlaps = paste(args.df[,1]+1,args.df[,2],sep="/")
  #Term Overlap   P.value Adjusted.P.value     Genes
  adjpvalues = p.adjust(pvalues, method = "BH")
  genes = lapply(intersects[rownames(args.df)],function(x){paste(x[[1]], collapse=";")})
  res <- data.frame(Term = terms,
                    Overlap = overlaps,
                    P.value = pvalues,
                    Adjusted.P.value = adjpvalues,
                    Genes = unlist(genes),
                    row.names = NULL,
                    stringsAsFactors = FALSE)
  return(res)
}


# human_line_plots produces plots based on 

human_line_plots = function(moduleDF, modulecol, regcol,
                     logFCs,  pwys, 
                     maxterms = 3,
                     clustercolors = NULL,
                     linetheme = theme(),
                     verbose = TRUE){
  
  if(is.null(clustercolors)){
    clustercolors = gg_color_hue(
      length(unique(as.character(moduleDF[,modulecol]))))
    names(clustercolors) = sort(unique(as.character(moduleDF[,modulecol])))
  }
  
  alllists = lapply(sort(unique(as.character(moduleDF[,modulecol]))), 
                    function(x){list()})
  names(alllists) = sort(unique(as.character(moduleDF[,modulecol])))
  
  alllists <- lapply(names(alllists),  function(i){
    mylist = lapply(c("up","down"), function(j){
      if(verbose){message("new KEGG")}
      moduleDF_module = moduleDF[moduleDF[,modulecol] %in% i,]
      listgenes =
        moduleDF_module[moduleDF_module[,regcol] == j,"mouse"]
      keggs = ownORA(listgenes, pathways = pwys,
                     universe = universe)
      keggs_sig = keggs[keggs$Adjusted.P.value<0.05,]
      keggs_sig = keggs_sig[order(keggs_sig$Adjusted.P.value),]
      keggs_sig$Term = gsub(" - Mus musculus \\(house mouse\\)","",keggs_sig$Term )
      
      listkeggs = keggs_sig
      
      xtype = substr(colnames(logFCs)[1],7,7)
      xbreaks = list(W = 0:7, D = seq(0,60,10))[[xtype]]
      xtitle = c(W = "Week", D = "Day")[xtype]
      
      logFCs = logFCs[listgenes,]
      logFCs$gene = rownames(logFCs)
      listplot = ggplot(reshape2::melt(logFCs), 
                        aes(y = value, 
                            x = as.numeric(gsub("logFC\\.[DW]","",variable)))) +
        geom_line(aes(group = gene), alpha = 0.4) +
        linetheme +
        geom_smooth(color = clustercolors[i], alpha = 0, method = "loess") +
        labs(x = xtitle, y = "Expression") +
        scale_x_continuous(breaks = xbreaks)
      return(list(listplot = listplot, listkeggs = listkeggs))
    })
    names(mylist )= c("up","down")
    return(mylist)
  })
  names(alllists) = sort(unique(as.character(moduleDF[,modulecol])))
  plotlist = lapply(alllists,
                    function(x){
                      lapply(x, "[[","listplot")})
  kegglist = lapply(alllists,
                    function(x){
                      lapply(x, "[[","listkeggs")})
  
  modules = names(plotlist)
  for(i in modules){
    print(plot_grid(ggdraw(textGrob(paste("Module",i), rot = 90)),
                    plot_grid(ggdraw(textGrob(
                      paste0("Up-regulated (",
                             length(unique(plotlist[[i]][["up"]]$data$gene)),")"),
                      rot = 90,
                      gp = gpar(fontface = "italic", fontsize = 8))),
                      ggdraw(textGrob(
                        paste0("Down-regulated (",
                               length(unique(plotlist[[i]][["down"]]$data$gene)),")"),
                        rot = 90, gp = gpar(fontface = "italic", fontsize = 8))),
                      nrow = 2),
                    plot_grid(plot_grid(plotlist = plotlist[[i]], nrow = 2),
                              plot_grid(plotlist = lapply(kegglist[[i]],function(j){
                                if(nrow(j) == 0){
                                  return(ggdraw(textGrob("No significant results",
                                                         gp = gpar(fontface = "italic", 
                                                                   fontsize = 8))))
                                }else{
                                  extra = 0
                                  if(nrow(j)>maxterms){extra = nrow(j)-maxterms;j = j[1:maxterms,]}
                                  j$Term = factor(j$Term, levels = rev(j$Term))
                                  
                                  topterms = paste(unlist(lapply(1:nrow(j), function(r){
                                    paste0(j$Term[r],
                                           " (",j$Overlap[r],")\n  ",
                                           paste(unlist(strsplit(j$Genes[r],";")),
                                                 collapse = " "))
                                  })),collapse="\n")
                                  if(extra>0){
                                    topterms = paste0(topterms, "\n+ ",extra," pathways")
                                  }
                                  return(ggdraw(textGrob(topterms,
                                                         just = "left", x = 0,
                                                         gp = gpar(col = "#555555", fontsize = 8)),
                                                clip = "on"))
                                }}), nrow = 2),
                              rel_widths = c(1,3), ncol = 2), ncol = 3, 
                    rel_widths = c(1,1,11)))
  }
}


combineandfilter = function(datalist){
  mergedcounts = datalist[[1]]
  if(length(datalist)>1){
    for(i in 2:length(datalist)){
      mergedcounts = merge(mergedcounts, datalist[[i]], by="row.names")
      rownames(mergedcounts) = mergedcounts$Row.names
      mergedcounts = mergedcounts[,-1]
    }
  }
  sel <- rowSums( mergedcounts > 5 ) >= ncol(mergedcounts)*0.15
  filt_counts = mergedcounts[sel,]
  filt_counts = filt_counts[,colSums(filt_counts)>0]
  return(filt_counts)
}

combinemeta = function(metalist,commoncns,commondata){
  if(length(metalist)>1){
    commonmeta = do.call("rbind",lapply(metalist,function(x){x[,commoncns]}))
  }else{
    commonmeta = metalist[[1]]
  }
  
  rownames(commonmeta) = commonmeta$SampleID
  commonmeta = commonmeta[match(colnames(commondata),rownames(commonmeta)),]
  
  commonmeta$nFeatures <- colSums(commondata>0)
  
  commonmeta$nCounts <- colSums(commondata)
  
  is_mito <- grepl("^MT-",rownames(commondata))
  commonmeta$perc_mito <- colSums(commondata[is_mito,]) / colSums(commondata) * 100
  
  is_ribo <- grepl("^RP[SL]",rownames(commondata))
  commonmeta$perc_ribo <- colSums(commondata[is_ribo,]) / colSums(commondata) * 100
  
  return(commonmeta)
}







getHumanDEs <- function(TCT_human, homologs){
  
  library(openxlsx)
  library(data.table)
  library(edgeR)
  
  Taman_UCDEG = 
    read.xlsx("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6290885/bin/jjx139_suppl_supplementary_table_1.xlsx", 
              startRow = 2)
  
  Taman_UCDEG = Taman_UCDEG[Taman_UCDEG$padj<0.01,]
  
  Taman_UC_mouse = homologs[match(Taman_UCDEG$Gene.Symbol, homologs$human),]
  Taman_UC_mouse = Taman_UC_mouse[!is.na(Taman_UC_mouse$mouse),]
  
  Taman_UC_mouse$logFC = Taman_UCDEG[match(Taman_UC_mouse$human,Taman_UCDEG$Gene.Symbol),"log2FoldChange"]
  
  ### Counts
  
  con <- gzcon(url(
    paste("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111889&format=file&file=GSE111889%5Fhost%5Ftx%5Fcounts%2Etsv%2Egz",
          sep="")))
  txt <- readLines(con)
  GSE111889counts <- read.table(textConnection(txt), header = TRUE, row.names = 1)
  
  con <- gzcon(url(
    paste("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE111nnn/GSE111889/matrix/GSE111889_series_matrix.txt.gz",
          sep="")))
  txt <- readLines(con)
  dat <- read.csv2(textConnection(txt), row.names = 1)
  
  ii <- readLines(con)
  tmp = t(read.delim( textConnection(txt) ,skip = (1:length(ii))[ii==""],header = F))
  
  colnames(tmp) <- gsub("[!]","",tmp[1,])
  rownames(tmp) <- tmp[,"Sample_geo_accession"]
  tmp <- tmp[-1,]
  GSE111889meta = data.frame(tmp)
  
  ### Metadata
  
  commoncns = c("SampleName","SampleID","Disease","Site","Specification","Study")
  
  GSE111889meta$Disease = as.character(c("UC","CD","Control")[match(substr(GSE111889meta$Sample_characteristics_ch1.1,17,40),
                                                                    c("Ulcerative Colitis","Crohn's Disease","Non IBD"))])
  GSE111889meta$SampleID = GSE111889meta$Sample_title
  rownames(GSE111889meta) = GSE111889meta$SampleID
  GSE111889meta$SampleName = GSE111889meta$ID_REF
  GSE111889meta$Site = gsub("\\(.*\\)","Colon",substr(GSE111889meta$Sample_characteristics_ch1,18,40))
  GSE111889meta$Specification = tolower(substr(GSE111889meta$Sample_characteristics_ch1.2,9,40))
  GSE111889meta$Study = "GSE111889"
  
  filt_counts = combineandfilter(list(
    GSE111889counts[,rownames(GSE111889meta[c(grep("ectum",GSE111889meta$Site),
                                              grep("olon",GSE111889meta$Site)),])]))
  commoncpm = log2(edgeR::cpm(filt_counts)+1)
  
  commonmeta = combinemeta(list(GSE111889meta), commoncns, filt_counts)
  
  dgList <- DGEList(counts=filt_counts, genes=rownames(filt_counts)) #table of counts
  dgList <- calcNormFactors(dgList)
  sampleType <- factor(commonmeta$Disease, levels=c("Control","CD","UC"))
  designMat <- model.matrix(~sampleType)
  dgList <- estimateGLMCommonDisp(dgList, design=designMat) #qCML common dispersion or tagwise
  #dispersions
  
  fit <- glmQLFit(dgList, designMat)
  de_CDvsCtrl = glmQLFTest(fit, coef = 2)
  de_CDvsCtrl_res = topTags(de_CDvsCtrl, n=nrow(fit))
  
  de_UCvsCtrl = glmQLFTest(fit, coef = 3)
  de_UCvsCtrl_res = topTags(de_UCvsCtrl, n=nrow(fit))
  
  de_UCvsCtrl_res_sig = de_UCvsCtrl_res$table[de_UCvsCtrl_res$table$FDR<0.01 & abs(de_UCvsCtrl_res$table$logFC)>1,]
  
  Lloyd_UC_mouse = homologs[match(rownames(de_UCvsCtrl_res_sig), homologs$human),]
  Lloyd_UC_mouse = Lloyd_UC_mouse[!is.na(Lloyd_UC_mouse$mouse),]
  
  Lloyd_UC_mouse$logFC = de_UCvsCtrl_res_sig[match(Lloyd_UC_mouse$human,rownames(de_UCvsCtrl_res_sig)),"logFC"]
  
  
  de_CDvsCtrl_res_sig = de_CDvsCtrl_res$table[de_CDvsCtrl_res$table$FDR<0.01 & abs(de_CDvsCtrl_res$table$logFC)>1,]
  
  Lloyd_CD_mouse = homologs[match(rownames(de_CDvsCtrl_res_sig), homologs$human),]
  Lloyd_CD_mouse = Lloyd_CD_mouse[!is.na(Lloyd_CD_mouse$mouse),]
  
  Lloyd_CD_mouse$logFC = de_CDvsCtrl_res_sig[match(Lloyd_CD_mouse$human,rownames(de_CDvsCtrl_res_sig)),"logFC"]
  
  DE_DF = unique(rbind(Lloyd_CD_mouse[,c(1,2)],Lloyd_UC_mouse[,c(1,2)], 
                       Taman_UC_mouse[,c(1,2)], TCT_human[,c(1,2)]))
  
  DE_DF$TCTmodule = TCT_human[match(DE_DF$mouse, TCT_human$mouse),"cluster"]
  
  DE_DF$Taman_UC_logFC = Taman_UC_mouse[match(DE_DF$human, Taman_UC_mouse$human), "logFC"]
  DE_DF$Lloyd_UC_logFC = Lloyd_UC_mouse[match(DE_DF$human, Lloyd_UC_mouse$human),"logFC"]
  DE_DF$Lloyd_CD_logFC = Lloyd_CD_mouse[match(DE_DF$human, Lloyd_CD_mouse$human),"logFC"]
  
  DE_DF$Taman_UC_reg = c(NA,"down","up")[sapply(DE_DF$Taman_UC_logFC, function(x){
    if(!is.na(x)){ (x>0)+2 }else{ 1 }})]
  DE_DF$Lloyd_UC_reg = c(NA,"down","up")[sapply(DE_DF$Lloyd_UC_logFC, function(x){
    if(!is.na(x)){ (x>0)+2 }else{ 1 }})]
  DE_DF$Lloyd_CD_reg = c(NA,"down","up")[sapply(DE_DF$Lloyd_CD_logFC, function(x){
    if(!is.na(x)){ (x>0)+2 }else{ 1 }})]
  
  return(DE_DF)
}

####