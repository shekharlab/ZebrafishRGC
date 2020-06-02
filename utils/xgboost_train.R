# Training a multiclass classifier
# input
#' train_Data: a gene expression matrix (variable genes x cells)
#' train_labels = factor of labels
#' var.genes = features to use for learning
#' do.scale = whether we want to z-score the features
#' train.frac = fraction of cells in each ident we want to use for training
XGBoost_train = function(train_Data, train_labels = NULL,var.genes=NULL, do.scale=FALSE,scale.mean = NULL, scale.var = NULL,max.cells.per.ident = 400, train.frac = 0.6){
  library(xgboost)
  if (!is.factor(train_labels)){
    train_labels = factor(train_labels)
  }
  
  if (length(train_labels) != ncol(train_Data)) stop("Error: There must be as many training IDs as there are cells (columns)")
  
  train_labels = train_labels[colnames(train_Data)]
  
  # Subsetting and Scaling
  train_Data = as.matrix(train_Data[var.genes,])
  if (do.scale){
    if (is.null(scale.mean) & is.null(scale.var)){
      scale.mean = rowMeans(train_Data)
      scale.var = apply(train_Data, 1, sd)
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
    } else {
      train_Data = t(scale(t(train_Data), center=scale.mean, scale=scale.var))
    }
    
  }
  
  # Training set vs. validation set
  training.set = c(); validation.set=c()
  training.label = c(); validation.label=c();
  print(paste0("Using either ", train.frac*100, " percent cells or ", max.cells.per.ident, " cells per cluster for training, whichever is smaller"))
  
  
  for (cc in c(1:length(levels(train_labels)))){
    i = levels(train_labels)[cc]
    
    cells.in.clust = names(train_labels)[train_labels == i];
    n = min(c(max.cells.per.ident, round(length(cells.in.clust)*train.frac)))
    train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
    validation.temp = setdiff(cells.in.clust, train.temp)
    training.set = c(training.set,train.temp); validation.set=c(validation.set,validation.temp)
    training.label = c(training.label, rep(cc-1,length(train.temp))); validation.label = c(validation.label, rep(cc-1, length(validation.temp)));
    # Consider upsampling
  }
  
  train_matrix <- xgb.DMatrix(data = t(train_Data[,training.set]), label=training.label)
  validation_matrix <- xgb.DMatrix(data = t(train_Data[,validation.set]), label=validation.label)
  
  numberOfClasses <- length(unique(training.label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "eta" = 0.2,"max_depth"=6, subsample = 0.6)
  nround    <- 200 # number of XGBoost rounds
  print(1)
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out validation set
  validation_pred <- predict(bst_model, newdata = validation_matrix)
  validation_prediction <- matrix(validation_pred, nrow = numberOfClasses,
                                  ncol=length(validation_pred)/numberOfClasses)
  
  valid_predlabels=apply(validation_prediction,2,which.max)-1
  A = table(validation.label, valid_predlabels)
  colnames(A) = levels(train_labels); rownames(A) = levels(train_labels)
  plotConfusionMatrix(A, order="Row", xlab.use = "True", ylab.use = "Predicted")
  
  to.return = list()
  to.return$bst_model = bst_model
  to.return$scale_mean = scale.mean
  to.return$scale_var = scale.var
  to.return$test_mat = A
  return(to.return)
}


# Plotting the Confusion matrix

plotConfusionMatrix = function(X,row.scale=TRUE, col.scale=FALSE, col.low="blue", col.high="red", max.size=5, ylab.use="Known", xlab.use="Predicted", order=NULL, x.lab.rot=FALSE, plot.return=TRUE, max.perc=100, title.use = "Validation test"){
  
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
  
  library(ggplot2)
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low =col.low,   high = col.high, limits=c(0, 100 ))+scale_size(range = c(1, max.size), limits = c(0,max.perc))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic")) + ggtitle(title.use) 
  
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  }
  print(p)
  
  if (plot.return) return(p)
}