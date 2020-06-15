
#' Creates a dotplot that shows how genes (columns) are expressed in different clusters (rows)
#'
#' @param object A Seurat object
#' @param genes.use A vector of genes to plot
#' @param use.counts A boolean indicating whether to use raw counts or scaled data
#' @param ident.use Which ident to use for plotting clusters
#' @param norm.exp Value between 0 and 1 to normalize expression to 
#' @param min.perc Minimum expression filter
#' @param max.val.perc Maximum percantage value: any clusters with a higher percentage are set to this value
#' @param max.val.exp Maximum expression value: any clusters with a higher expression are set to this value
#' @param make.clust.diag Boolean of whether to diagonalize the gene order or not
#' @param do.plot Whether to construct a plot (TRUE) or just return the percentage and expression matrices (FALSE)
#' @param clusters_as_rows Whether to plot clusters as rows (TRUE) or columns (FALSE)
#' @param alpha.use Alpha value to use for ggplot
#' @param col.low Color for low expression values
#' @param col.high Color for high expression values
#' @param max.size Maximum dot size
#' @param title.use Main title
#'
#' @examples
#' basic_dotplot(pbmc, genes.use = c("IL32", "LTB", "CFD"), col.high = "blue", make.clust.diag = TRUE)
basic_dotplot = function(object,genes.use=NULL, use.counts=FALSE, ident.use = NULL, 
                         norm.exp=NULL, min.perc=0,
                         max.val.perc=NULL, max.val.exp=NULL, 
                         make.clust.diag = FALSE, do.plot=TRUE,
                         clusters_as_rows = TRUE, alpha.use = 1,
                         col.low="white", col.high="red",
                         max.size=10, title.use=NULL) {
  if (use.counts){
    data.use = object@assays$RNA@counts
    data.use = data.use[intersect(genes.use, rownames(data.use)),]
  } else {
    data.use = object@assays$RNA@data
    data.use = data.use[intersect(genes.use, rownames(data.use)),]
    data.use = exp(data.use) - 1
  }
  
  
  # Choose clusters
  if (is.null(ident.use)){
    ident.use = levels(Idents(object))
  } else {
    ident.use = intersect(ident.use, levels(Idents(object)))
  }
 
  # Prepare data for plotting
  PercMat = matrix(0, nrow=length(genes.use), ncol = length(ident.use))
  rownames(PercMat) = genes.use; colnames(PercMat) = ident.use
  
  #Matrix of average transcript levels
  ExpMat = PercMat;
  
  # Populate PercMat and ExpMat
  for (i in ident.use){
    cells_in_cluster = WhichCells(object, idents = i)
    PercMat[,i] = tfPercentExpression(object, i, genes.use, threshold = 0)
    vec.exp = apply(data.use[, cells_in_cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {0})
    ExpMat[,i] = vec.exp
  }
  
  # Normalize expression
  if (!is.null(norm.exp)){
    if (norm.exp < 0 | norm.exp > 1){
      print("Warning: norm.exp should be a value between (0,1). Skipping normalization")
      next
    } else{
      quant.vals = apply(ExpMat,1, function(x) quantile(x, norm.exp))
      ExpMat = t(scale(t(ExpMat), center=FALSE, scale=quant.vals))
    }
  }
  
  # Apply minimum expression filter
  genes.use_filt = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
  PercMat = PercMat[genes.use_filt,]
  ExpMat = ExpMat[genes.use_filt,]
  
  # Further thresholding of perc and exp
  if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
  if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
  
  # Code to make output diagonal
  if (make.clust.diag){
    library(gplots)
    hc <- heatmap.2(PercMat, hclustfun = function(x) hclust(x,method="single"))
    PercMat <- PercMat[hc$rowInd, hc$colInd]
    ExpMat <- ExpMat[hc$rowInd, hc$colInd]
  }  
  
  
  if (do.plot){
    
    ExpVal = melt(ExpMat)
    PercVal = melt(PercMat)
    colnames(ExpVal) = c("gene","cluster","Exp")
    ExpVal$percExp = PercVal$value*100
    ExpVal$gene = factor(ExpVal$gene, levels=genes.use)
    ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
    
    if (clusters_as_rows){
      # Initiate dotplot
      p=ggplot(ExpVal, aes(y = cluster,  x = gene)) + geom_point(aes(colour = Exp,  size =percExp), alpha=alpha.use)  
      # Add colors
       p = p + scale_color_gradient(low = col.low,   high = col.high, limits=c( 0, max(ExpVal$Exp) )) + scale_size(range = c(0, max.size))+   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
       
       # Add title and labels
       p = p + ggtitle(title.use) + theme(plot.title = element_text(hjust = 0.5, size=15))
      p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) +  theme(axis.text.y=element_text(size=12, face="italic"))
    } else {
      ExpVal$gene = factor(ExpVal$gene, levels=genes.use)
      ExpVal$cluster = factor(ExpVal$cluster, levels= rev(ident.use))
      
      # Initiate dotplot
      p=ggplot(ExpVal, aes(y = gene,  x = cluster)) + geom_point(aes(colour = Exp,  size =percExp), alpha=alpha.use)  
      # Add colors
      p = p + scale_color_gradient(low = col.low,   high = col.high, limits=c( 0, max(ExpVal$Exp) )) + scale_size(range = c(0, max.size))+   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      # Add title and labels
      p = p + ggtitle(title.use) + theme(plot.title = element_text(hjust = 0.5, size=15))
      p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) +  theme(axis.text.y=element_text(size=12, face="italic"))
    }
    
    print(p)
    
  } else {
    return(list(PercMat, ExpMat))
  }
  

}


#' Plots a confusion matrix
#'
#' @param X A confusion matrix
#' @param row.scale If TRUE, scales the rows to sum to one
#' @param col.scale If TRUE, scales the columns to sum to one. Both row.scale and col.scale cannot be TRUE
#' @param col.low Color to use for low values
#' @param col.high Color to use for high values
#' @param max.size Maximum size of dots
#' @param ylab.use Y axis label
#' @param xlab.use X axis label
#' @param order Set to either "Row" or "Col" to order the plot according to the row idents or column idents. By default, it will plot in order of the rownames of the Confusion matrix
#' @param x.lab.rot If TRUE, rotate the labels of the x-axis
#' @param plot.return Whether to return the plot or not
#' @param max.perc Maximum percentage. Any dot with a percentage higher will be set to this value. 
#' 
#' @return A ggplot object of the confusion matrix
plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=FALSE, plot.return=TRUE, max.perc=100){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  X[is.na(X)] = 0
  if (max(X) > 100){
    X=X/100
  }
  
  orig.rownames = rownames(X)
  orig.colnames = colnames(X)
  if (!is.null(order)){
    if (order == "Row"){  
      factor.levels = c()
      for (i1 in colnames(X)){
        if (max(X[,i1]) < 50) next
        ind.sort = rownames(X)[order(X[,i1], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
    
    if (order == "Col") {
      factor.levels = c()
      for (i1 in rownames(X)){
        if (max(X[i1,]) < 50) next
        ind.sort = rownames(X)[order(X[i1,], decreasing=TRUE)]
        ind.sort = ind.sort[!(ind.sort %in% factor.levels)]
        factor.levels = c(factor.levels, ind.sort[1])
      }
      factor.levels = c(factor.levels, setdiff(rownames(t(X)), factor.levels))
      factor.levels = factor.levels[!is.na(factor.levels)]
    } 
  } else {
    factor.levels = rownames(t(X))
  }
  
  factor.levels = c(factor.levels, setdiff(rownames(X), factor.levels))
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  #X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  #X$Predicted = factor(X$Predicted, levels = rev(factor.levels))
  
  if (!is.null(order)){
    if (order == "Row"){ 
      X$Known = factor(X$Known, levels=rev(factor.levels));
      X$Predicted = factor(X$Predicted, levels = orig.colnames)
      
    }
    if (order == "Col"){
      X$Predicted = factor(X$Predicted, levels = factor.levels);
      X$Known = factor(X$Known, levels=rev(orig.rownames));
    }
  } else {
    X$Known = factor(X$Known, levels=rev(unique(X$Known)));
    X$Predicted = factor(X$Predicted, levels=unique(X$Predicted));
  }
  
  #print(sum(is.na(X$Known)))
  
  
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =col.low,   high = col.high, limits=c(0, 100 ))+scale_size(range = c(1, max.size), limits = c(0,max.perc))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))  
  
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  print(p)
  
  if (plot.return) return(p)
}


#' Generates a diagonlized confusion matrix plot
#'
#' @param C A Confusion matrix
#' @param xlab.use X label to use
#' @param ylab.use Y label to use
#' @param title.use Main title
#' @param col.low Color to use for low values
#' @param col.high Color to use for high values
#' @param order Whether to return the order of clusters plotted on the x axis
#' @return A confusion matrix plot as well as the row order
MakePrettyConfusionMatrix = function(C, xlab.use = "Training clusters", ylab.use = "Test clusters",
                                     title.use = "Performance Confusion Matrix", col.high = "darkblue",
                                     col.low = "white", order = FALSE){
  library(gplots)
  hc <- heatmap.2(C,hclustfun = function(x) hclust(x,method="single"))
  C1=C[hc$rowInd, hc$colInd]
  row.max = apply(C1,1,which.max)
  names(row.max) = rownames(C1)
  row.ord = names(sort(row.max))
  
  plotConfusionMatrix(C1[row.ord,], xlab.use = xlab.use, ylab.use = ylab.use, col.high = col.high, col.low = col.low)
  
  if(order){
    return(row.ord)
  }
}


#' Plots the relative cluster size of each cluster in an object
#'
#' @param A Seurat object
PlotClusterSize = function(object){
  num_cells <- dim(object@meta.data)[1]
  cluster_percent <- vector()
  for(i in levels(object@meta.data$clusterID)){
    clus_cells <- length(which(object@meta.data$clusterID == i))
    cluster_percent[i] <- clus_cells / num_cells * 100
  }
  barplot(cluster_percent, xlab = "Cluster ID", ylab = "Percentage of Cells [%]")
}

#' Given a dataframe generated by FindAllMarkers, plots the top DE gene for each cluster and returns
#' the list of genes plotted. If edits = TRUE, the user can manually choose which gene is plotted 
#' for which cluster by selecting from the top 10 differentially expressed genes.
#'
#' @param object A Seurat object
#' @param markers A dataframe of DE genes generated by FindAllMarkers
#' @param edits Boolean value indicating whether manual edits should be made to the plot
#' @param plot Boolean value indicating whether the plot should be shown or not.
#'
#' @return A vector of genes plotted
#' @examples
#' PlotUniqueMarkers(object = pbmc, markers = pbmc.markers)
#' PlotUniqueMarkers(object = pbmc, markers = pbmc.markers, edits = TRUE)
PlotUniqueMarkers = function(object, markers, edits = FALSE, plot = TRUE){
  plot_markers <- NULL
  for(i in levels(Idents(object))){
    index <- 1
    top_mark <- head(subset(markers, cluster == i))[index,]$gene
    while(top_mark %in% plot_markers){
      index <- index+1
      top_mark <- head(subset(markers, cluster==i))[index,]$gene
    }
    plot_markers <- c(plot_markers, top_mark)
  }
  
  if (plot){
    plot <- DotPlot(object, features = plot_markers) + RotatedAxis()
    print(plot)
  } 
  if (edits){
    clust <- readline(prompt = "Change which cluster? Press q to exit: ")
    while(clust != "q"){
      if(clust %in% levels(Idents(object))){
        indent <- match(clust, levels(Idents(object)))
      }
      else{
        break()
      }
      subplot <- DotPlot(object, features = head(subset(markers, cluster == clust)$gene, 10)) + RotatedAxis()
      print(subplot)
      new_gene <- readline(prompt = "Use which gene? ")
      plot_markers[indent] <- new_gene
      plot <- DotPlot(object, features = plot_markers) + RotatedAxis()
      print(plot)
      clust <- readline(prompt = "Change which cluster? Press q to exit: ")
    }
  }
  return(plot_markers)
}

