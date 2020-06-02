
#' Converts gene names within a Seurat object to lowercase
#'
#' @param object A Seurat object.
#' @param integration Boolean value indicating whether the object has an integrated assay.
#'
#' @return Returns a Seurat object with genes in lowercase.
LowerCase_genes = function(object, integration = FALSE){
  
  rownames(object@assays$RNA@counts)<- tolower(rownames(object@assays$RNA@counts))
  rownames(object@assays$RNA@data)<- tolower(rownames(object@assays$RNA@data))
  rownames(object@assays$RNA@scale.data)<- tolower(rownames(object@assays$RNA@scale.data))
  object@assays$RNA@var.features <- tolower(object@assays$RNA@var.features)

  if (integration){
    rownames(object@assays$integrated@counts)<- tolower(rownames(object@assays$integrated@counts))
    rownames(object@assays$integrated@data)<- tolower(rownames(object@assays$integrated@data))
    rownames(object@assays$integrated@scale.data)<- tolower(rownames(object@assays$integrated@scale.data))
    object@assays$integrated@var.features <- tolower(object@assays$integrated@var.features)
  }
  
  return(object)
}


#' Converts gene names withing a Seurat object to uppercase 
#'
#' @param object A Seurat object.
#' @param integration Boolean value indicating whether the object has an integrated assay.
#'
#' @return Returns a Seurat object with genes in uppercase.
UpperCase_genes = function(object, integration = FALSE){
  
  rownames(object@assays$RNA@counts)<- toupper(rownames(object@assays$RNA@counts))
  rownames(object@assays$RNA@data)<- toupper(rownames(object@assays$RNA@data))
  rownames(object@assays$RNA@scale.data)<- toupper(rownames(object@assays$RNA@scale.data))
  object@assays$RNA@var.features <- toupper(object@assays$RNA@var.features)
  #rownames(object@assays$RNA@meta.features)<- toupper(rownames(object@assays$RNA@meta.features))
  
  if (integration){
    rownames(object@assays$integrated@counts)<- toupper(rownames(object@assays$integrated@counts))
    rownames(object@assays$integrated@data)<- toupper(rownames(object@assays$integrated@data))
    rownames(object@assays$integrated@scale.data)<- toupper(rownames(object@assays$integrated@scale.data))
    object@assays$integrated@var.features <- toupper(object@assays$integrated@var.features)
    #rownames(object@assays$integration@meta.features)<- toupper(rownames(object@assays$integration@meta.features))
  }
  
  return(object)
}


#' This function calculates the percentage of cells within a cluster that express a certain gene
#' 
#' @param A Seurat object
#' @param clusID Cluster of interest
#' @param tf Gene of interest
#' @param threshold Threshold above which a gene must be expressed to be counted as expressed within that cell (default 0)
#'
#' @return The percentage of cells within a cluster that express the gene
#'
#' @examples
#' tfPercentExpression(object = pbmc, clusID = 1, tf = "LTB")
tfPercentExpression = function(object, clusID, tf, threshold = 0){
  # Determine which cells are in each cluster
  cluster_cells <- WhichCells(object, idents = clusID)
  
  # Determine which cells express transcription factors over a given threshold
  positiveCellMat <- object@assays$RNA@data[tf,cluster_cells] > threshold
  
  # Calculate the percentage of cells that express the transcription factors
  if(class(positiveCellMat) == "logical"){
    numPositiveCells <- sum(positiveCellMat)
  }else{
  numPositiveCells <- Matrix::rowSums(positiveCellMat)
  }
  numCells <- length(cluster_cells)
  percent_express <- numPositiveCells / numCells
  return(percent_express)
}


#' Calculate the Shannon Entropy for a vector
#'
#' @param vector 
#' @param threshold A threshold set to remove background. When the vector is scaled to sum to one, any elements below the threshold will be set to zero and the vector will be scaled.
#'
#' @return The Shannon Entropy for the vector
#'
#' @examples
#' ShannonEntropy(vector = c(1, 0, 0, 1, 0))
#' ShannonEntropy(vector = c(10, .5, 0, 10, 0), threshold = .05)
ShannonEntropy = function(vector, threshold = 0){
  # Scale the vector by removing background
  vector <- vector / sum(vector)
  vector[vector<threshold] <- 0
  vector <- vector / sum(vector)
  # Calculate Shannon entropy
  Hx <- -sum(vector[vector>0]*log2(vector[vector>0]))
  return(2^Hx)
}

#' Calculate the Occupation Number for a vector
#'
#' @param vector 
#' @param threshold A threshold set to remove background. When the vector is scaled to sum to one, any elements below the threshold will be set to zero and the vector will be scaled.
#'
#' @return The Occupation Number for the vector
#'
#' @examples
#' OccupationNumber(vector = c(1, 0, 0, 1, 0))
#' OccupationNumber(vector = c(10, .5, 0, 10, 0), threshold = .05)
OccupationNumber = function(vector, threshold = 0){
  # Scale the vector by removing background
  vector <- vector / sum(vector)
  vector[vector<threshold] <- 0
  vector <- vector / sum(vector)
  
  return(1/sum(vector^2))
  
}


#' Finds genes that are uniquely expressed in a single cluster from a list of differentially expressed genes
#'
#' @param object A Seurat object
#' @param markers A dataframe generated by FindAllMarkers
#' @param DEgenes_to_check Number of DE genes to check for uniqueness
#' @param min_percent_expression If a gene is not expressed in a cluster above this threshold, it is not tested for uniqueness
#'
#' @return A vector of genes that are uniquely expressed in each cluster
#'
#' @examples
#' FindUniqueMarkers(object = pbmc, markers = pbmc.markers)
FindUniqueMarkers = function(object, markers, DEgenes_to_check = 30, min_percent_expression = .3){
  object <- UpperCase_genes(object)
  markers$gene <- toupper(markers$gene)
  unique_markers <- NULL
  
  for(i in levels(Idents(object))){
    low_enrichment_score <- Inf
    top_mark <- NULL
    top_genes <- head(subset(markers, cluster == i), DEgenes_to_check)$gene
    
    for(gene in top_genes){
      per_express_i <- tfPercentExpression(object, clusID = i, tf = gene, threshold = 0)
      if(per_express_i < min_percent_expression){
        next()
      }
      per_express_vec <- unlist(lapply(levels(Idents(object)), function(x) tfPercentExpression(object, clusID = x, tf=gene, threshold = 1)))
      enrichment_score <- sum(per_express_vec / per_express_i)
      if(enrichment_score < low_enrichment_score){
        low_enrichment_score = enrichment_score
        top_mark = gene
        }
      }
    unique_markers <- c(unique_markers, top_mark)
    }
    
  return(unique_markers)
}

#' This function calculates the percentage of cells within all clusters that express a certain gene
#'
#' @param object A Seurat object
#' @param gene Gene (or genes) of interest
#' @param threshold Threshold above which a gene must be expressed to be counted as expressed within that cell (default 0)
#'
#' @return A matrix of genes by expression level in each cluster
#'
#' @examples
#' ExpressionByCluster(object = pbmc, gene = "LTB")
ExpressionByCluster = function(object, gene, threshold = 0){
  num_idents <- length(levels(Idents(object)))
  per_express <- lapply(1:num_idents, function(x) tfPercentExpression(object, clusID = x, tf=gene, threshold = threshold))
  per_express_mat <- matrix(unlist(per_express), ncol=num_idents, dimnames = list(gene, 1:num_idents))
  return(per_express_mat)
}


#' Simplified workflow for clustering a Seurat object
#'
#' @param object A Seurat object 
#' @param nfeatures Number of features to select as top variable features
#' @param numPCs Number of principal components for clustering analysis
#' @param normalization.method Method for normalization (select either LogNormalize, CLR, or RC). See Seurat::NormalizeData for more information.
#' @param scale.factor Scale factor for cell-level normalization
#' @param selection.method How to choose top variable features (select either vst, mean.var.plot, or dispersion). See Seurat::FindVariableFeatures for more information.
#' @param scale.all Boolean indicating whether to scale all features, or just variably expressed features
#' @param cluster_resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param elbow Boolean indicating whether to show an Elbow Plot for better selection of numPCs
#'
#' @return A clustered Seurat object with tSNE and UMAP embeddings
ClusterSeurat = function(object, nfeatures = 2000, numPCs = 30, normalization.method = "LogNormalize", scale.factor = 10000, selection.method = "vst",  scale.all = FALSE,  cluster_resolution = .5, elbow = FALSE){
  object <- NormalizeData(object, normalization.method = normalization.method, scale.factor = scale.factor)
  object <- FindVariableFeatures(object, selection.method = selection.method, nfeatures = nfeatures)
  if(scale.all){
    all.genes <- rownames(object)
    object <- ScaleData(object, features = all.genes)
  }
  else{
    object <- ScaleData(object)
  }

  object <- RunPCA(object)
  if(elbow){
    ElbowPlot(object, ndims = 50)
    numPCs <- readline(prompt = "Enter the number of PCs desired for clustering analysis: ")
  }
  object <- FindNeighbors(object, dims = 1:numPCs)
  object <- FindClusters(object, resolution = cluster_resolution)
  object <- RunTSNE(object, dims = 1:numPCs)
  object <- RunUMAP(object, dims = 1:numPCs)
  
  return(object)
}


#' Builds a confusion matrix given an xgboost classification model and a Seurat object to apply the model to 
#'
#' @param train A Seurat object that was used to train the xgboost model
#' @param test A Seurat object to apply the xgboost model to
#' @param model An xgboost model
#' @param scale.by.model If TRUE, the test dataset is scaled by the model mean and variance. If FALSE, the test dataset is z-scored.
#'
#' @return A Confusion matrix
#'
#' @examples
#' adult_to_larva <- BuildConfusionMatrix(larva, adult, model = larva.model)
BuildConfusionMatrix = function(train, test, model, scale.by.model = FALSE ){
  genes.use <- toupper(model$bst_model$feature_names)
  train <- UpperCase_genes(train)
  test <- UpperCase_genes(test)
  train_data = train@assays$RNA@data[genes.use,]
  train_id = Idents(train)
  test_data = test@assays$RNA@data[genes.use,]
  test_data = as.matrix(test_data)
  test_id = Idents(test)
  
  # Use trained model to predict on test data
  numberOfClasses = length(levels(train_id))
  test_xgb = t(test_data[genes.use,])
  if (scale.by.model){
    test_xgb = scale(test_xgb, center=model$scale_mean, scale = model$scale_var)
  }
  else {
    test_xgb = scale(test_xgb)
  }
  test_xgb = xgboost::xgb.DMatrix(test_xgb)
  
  test_pred <- predict(model$bst_model, newdata = test_xgb)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses)
  
  # Find best class for each cell
  test_pred_margins = apply(test_prediction,2,max)
  test_predlabels = apply(test_prediction,2,which.max)
  names(test_predlabels) = colnames(test_data)
  test_pred_names = levels(train_id)[test_predlabels]
  names(test_pred_names) = names(test_predlabels)
  
  confusion_matrix = table(test_id, test_pred_names)
  return(confusion_matrix)
}


#' Trains a xgboost classification model
#'
#' @param object A Seurat object to learn
#' @param training_genes A vector of training genes
#' @param train_ident Which ident for the model to learn. 
#' @param do.scale Boolean value indicating whether to scale the data or not.
#'
#' @return An xgboost model
TrainModel = function(object, training_genes, train_ident = NULL, do.scale = TRUE){
  object <- UpperCase_genes(object)
  training_genes <- toupper(training_genes)
  train_data = object@assays$RNA@data[training_genes,]
  if(is.null(train_ident)){
    train_id = Idents(object)
  }
  else{
    Idents(object) <- train_ident
    train_id <- Idents(object)
  }
  bst_model = XGBoost_train(train_data, train_labels = train_id, 
                                  var.genes = training_genes, 
                                  do.scale = do.scale)
  return(bst_model)
}

#' Diagonlizes a set of genes for a Dot Plot
#'
#' @param genes A vector of genes to plot
#' @param object A Seurat object
#' @param increasing Boolean indicating the direction of the diagonal
#'
#' @return The same vector of genes but in diagonal order
DiagonalizeGenes <- function(genes, object, increasing = FALSE){
  num_idents <- length(levels(Idents(object)))
  per_express <- lapply(1:num_idents, function(x) tfPercentExpression(object, clusID = x, tf=genes, threshold = 0))
  per_express_mat <- matrix(unlist(per_express), ncol=num_idents, dimnames = list(genes, 1:num_idents))
  per_express_mat <- per_express_mat[,as.numeric(levels(Idents(object)))]
  gene_max_cluster <- apply(per_express_mat, 1, which.max)
  order <- names(sort(gene_max_cluster))
  if(increasing){
    order <- rev(order)
  }
  return(order)
}

#' Constructs a dendrogram and creates a new column in the metadata that is factored according to the dendrogram
#'
#' @param object A Seurat object with a new column "dendro_order" that is factored by a dendrogram
DendroOrder <- function(object){
  tree_obj = object@tools$BuildClusterTree
  left_clusts = Seurat:::GetLeftDescendants(tree_obj, length(levels(object@meta.data$clusterID))+1)
  right_clusts = Seurat:::GetRightDescendants(tree_obj, length(levels(object@meta.data$clusterID))+1)
  clust_order = c(left_clusts, right_clusts)
  object@meta.data$dendro_order = factor(object@meta.data$clusterID, levels = clust_order)
  return(object)
}