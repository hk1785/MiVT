BoRT <-
function(Z, tax.prop, tax.rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), 
                 minsplit = 10, minbucket = 5, cp = 0.01, tree.pruning = FALSE, n.boot = 20000) {
  
  Z <- Z*100
  tax <- qc.out$tax.prop[tax.rank][[1]]
  
  XX <- as.data.frame(tax*100)
  taxa <- colnames(XX)
  colnames(XX) <- paste(substr(tax.rank, 1, 1), 1:ncol(XX), sep = "")
  ZX <- cbind(Z, XX)
  colnames(ZX)[1] <- "Y"
  
  par(mfrow = c(1,1))
  rpart.para <- rpart.control(minsplit = minsplit, minbucket = minbucket, cp = cp)
  tree.fit <- rpart(Y~., data = ZX, control = rpart.para)
  
  if (tree.pruning) {
    
    tree.cont <- tree::tree.control(nobs = nrow(ZX), mincut = minbucket, minsize = minsplit, mindev = cp)
    fit <- tree::tree(Y ~ ., data = ZX, method = "recursive.partition", split = "deviance", control = tree.cont)
    cv.fit <- tree::cv.tree(object = fit, K = n.folds)
    (best.num.leaves <- cv.fit$size[which.min(cv.fit$dev)])
    
    if (best.num.leaves == 0) {
      
      stop("no subgroups are found.")
      
    } else {
      
      num.leaves <- numeric()
      for (i in 1:100) {
        rpart.para <- rpart.control(minsplit = minsplit, minbucket = minbucket, cp = cp - 0.005 + 0.005*i)
        tree.fit.tri <- rpart(Y~., data = ZX, control = rpart.para)
        num.leaves[i] <- sum(tree.fit.tri$frame$var == "<leaf>")
      }
      ind <- which(num.leaves == best.num.leaves)[1]
      rpart.para <- rpart.control(minsplit = minsplit, minbucket = minbucket, cp = cp - 0.005 + 0.005*ind)
      tree.fit <- rpart(Y~., data = ZX, control = rpart.para)
      
    }
    
  }
  
  rpart.plot(tree.fit, type = 4, digits = 3, clip.right.labs = FALSE, extra = 0)
  ind <- grep(substr(tax.rank, 1, 1), as.character(tree.fit$frame$var))
  sel.taxa <- as.character(tree.fit$frame$var)[ind]
  taxa.num <- as.numeric(gsub(substr(tax.rank, 1, 1), "", as.character(tree.fit$frame$var)[ind]))
  Sel.Taxa <- t(as.data.frame(taxa[taxa.num]))
  rownames(Sel.Taxa) <- "Full names"
  colnames(Sel.Taxa) <- sel.taxa
  
  BoRT.out <- boot.test(Z.hat = Z, tree.fit, n.boot = n.boot)
  BoRT.out <- round(BoRT.out, 3)
  colnames(BoRT.out) <- NULL
  
  out <- list(Sel.Taxa = Sel.Taxa, BoRT.out = BoRT.out)
  
  return(out)
}
