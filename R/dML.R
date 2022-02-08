dML <-
function(y, Tr, X, tree, n.folds = 10, n.rep = 2, 
                      alpha = seq(0.05, 0.95, 0.05), n.trees = 1000, n.neus = c(1/2, 1/3, 1/4), 
                      interac = TRUE, multi.rarefy = 1, n.seed = 123) {
  
  set.seed(n.seed)
  
  X <- as.data.frame(X)
  X <- t(X)
  n <- nrow(X)
  
  if (multi.rarefy == 1) {
    n <- nrow(X)
    X.ra <- rarefy_even_depth(otu_table(X, taxa_are_rows = FALSE), min(rowSums(X)), rngseed = n.seed)
  } else {
    X.ra.list <- list()
    for (i in 1:multi.rarefy) {
      X.ra.list[[i]] <- rarefy_even_depth(otu_table(X, taxa_are_rows = FALSE), min(rowSums(X)), rngseed = i)
    }
    
    rare.taxa <- colnames(X.ra.list[[1]])
    for (i in 2:multi.rarefy) {
      rare.taxa <- intersect(rare.taxa, colnames(X.ra.list[[i]]))
    }
    
    ind <- match(rare.taxa, colnames(X.ra.list[[1]]))
    sum.X.ra<- X.ra.list[[1]][,ind]
    for (i in 2:multi.rarefy) {
      ind <- match(rare.taxa, colnames(X.ra.list[[i]]))
      sum.X.ra <- sum.X.ra + X.ra.list[[i]][,ind]
    }
    X.ra <- round(sum.X.ra/multi.rarefy)
  }
  
  yTs <- paste(y, Tr, sep = ".")
  len.yTs <- table(yTs)
  uni.yTs <- names(len.yTs)
  fol.len.yTs <- ceiling(len.yTs/n.folds)
  
  ind.folds <- create_folds(yTs, k = n.folds, m_rep = n.rep, type = "stratified")
  
  prop <- X/rowSums(X)
  imp.prop <- prop + 0.5
  X.clr <- compositions::clr(imp.prop)
  
  Dis.list <- list()
  D.uni <- GUniFrac(X, tree, alpha = c(0.5, 1))$unifracs
  Dis.list[[1]] <- as.matrix(proxy::dist(X.clr, method = "Euclidean"))
  Dis.list[[2]] <- as.matrix(proxy::dist(X.ra, method = "Jaccard"))
  Dis.list[[3]] <- as.matrix(bcdist(X.ra))
  Dis.list[[4]] <- D.uni[, , "d_UW"]
  Dis.list[[5]] <- D.uni[, , "d_0.5"]
  Dis.list[[6]] <- D.uni[, , "d_1"]
  
  X.coor <- list()
  
  for (q in 1:length(Dis.list)) {
    
    D <- Dis.list[[q]]
    C <- diag(n) - (1/n)*matrix(rep(1, n*n), n)
    B <- (-1/2) * C %*% D^2 %*% C
    ei <- eigen(B)
    ind <- which(ei$values > 0)
    X.coor[[q]] <- as.data.frame(ei$vectors[,ind] %*% diag(sqrt(ei$values[ind])))
    
  }
  
  dat.list <- list()
  
  for (q in 1:length(Dis.list)) {
    
    X <- X.coor[[q]]
    
    if (interac) {
      
      T0X <- abs(Tr-1)*X
      colnames(T0X) <- paste("T0", colnames(X), sep = "")
      
      T1X <- Tr*X
      colnames(T1X) <- paste("T1", colnames(X), sep = "")
      
      Tr.X <- cbind(Tr, X, T0X, T1X)
      colnames(Tr.X)[1] <- "T" 
      
      dat <- as.data.frame(cbind(y, Tr.X))
      
      Tr.X.1 <- cbind(rep(1, nrow(X)), X, rep(0, nrow(X))*X, rep(1, nrow(X))*X)
      colnames(Tr.X.1)[1] <- "T" 
      colnames(Tr.X.1)[2:ncol(Tr.X.1)] <- c(colnames(X), paste("T0", colnames(X), sep = ""), paste("T1", colnames(X), sep = ""))
      
      Tr.X.0 <- cbind(rep(0, nrow(X)), X, rep(1, nrow(X))*X, rep(0, nrow(X))*X)
      colnames(Tr.X.0)[1] <- "T" 
      colnames(Tr.X.0)[2:ncol(Tr.X.0)] <- c(colnames(X), paste("T0", colnames(X), sep = ""), paste("T1", colnames(X), sep = ""))
      
    } else {
      
      dat <- as.data.frame(cbind(y, Tr, X))
      Tr.X.1 <- cbind(rep(1, nrow(X)), X)
      colnames(Tr.X.1)[1] <- "T" 
      Tr.X.0 <- cbind(rep(0, nrow(X)), X)
      colnames(Tr.X.0)[1] <- "T" 
      
      y <- dat[,1]
      Tr.X <- dat[,-1]
      colnames(Tr.X)[1] <- "T"
      
    }
    
    dat.list[[q]] <- dat
    
  }
  
  nn.cv.err.dis <- list()
  nn.cv.cro.dis <- list()   
  
  rf.cv.err.dis <- list()
  rf.cv.cro.dis <- list()
  
  en.cv.err.dis <- list()
  en.cv.cro.dis <- list()
  
  en.para.dis <- list()
  
  for (q in 1:length(Dis.list)) {
    
    nn.cv.errs <- list()
    nn.cv.cros <- list()
    
    rf.cv.errs <- list()
    rf.cv.cros <- list()
    
    en.cv.errs <- list()
    en.cv.cros <- list()
    
    en.paras <- list()
    
    for (i in 1:(n.folds*n.rep)) {
      
      dat <- dat.list[[q]]
      dat.trai <- dat[ind.folds[[i]],]
      dat.test <- dat[-ind.folds[[i]],]
      
      # NN
      
      if (interac) {
        n.pred <- (ncol(dat.trai)-2)/3
      }
      n.neu <- ceiling(n.pred*n.neus)
      
      nn.cv.err <- numeric()
      nn.cv.cro <- numeric()
      
      for (j in 1:length(n.neu)) {
        
        nn.fit <- neuralnet(y ~ ., data = dat.trai, hidden = c(n.neu[j], ceiling(n.neu[j]/2), ceiling(n.neu[j]/4)), 
                            err.fct = "ce", act.fct = "logistic", algorithm = "rprop+", linear.output = FALSE)
        
        pred.val <- as.numeric(neuralnet::compute(nn.fit, dat.test[,-1])$net.result)
        
        com <- dat.test[,1] - (pred.val > 0.5)
        pred.val <- round(pred.val, 3)
        ind.1 <- which(pred.val == 1)
        ind.0 <- which(pred.val == 0)
        pred.val[ind.1] <- 0.999
        pred.val[ind.0] <- 0.001
        nn.cv.err[j] <- sum(com != 0)/nrow(dat.test)
        nn.cv.cro[j] <- - (1/nrow(dat.test)) * sum(dat.test[,1]*log(as.numeric(pred.val)) + (1-dat.test[,1])*log(1-as.numeric(pred.val)))
        
      }
      
      nn.cv.errs[[i]] <- nn.cv.err
      nn.cv.cros[[i]] <- nn.cv.cro   
      
      # RF
      
      if (interac) {
        n.pred <- (ncol(dat.trai)-2)/3
      }
      mtry <- ceiling(c(n.pred*(1/2), sqrt(n.pred), log(n.pred)))
      
      rf.cv.err <- numeric()
      rf.cv.cro <- numeric()
      
      for (j in 1:length(mtry)) {
        
        rf.fit <- randomForest(y = as.factor(dat.trai$y), x = dat.trai[,-1], mtry = mtry[j], ntree = n.trees, family = "binomial")
        
        pred.val <- as.numeric(predict(rf.fit, dat.test[,-1], "prob")[,2])
        com <- dat.test[,1] - (pred.val > 0.5)
        pred.val <- round(pred.val, 3)
        ind.1 <- which(pred.val == 1)
        ind.0 <- which(pred.val == 0)
        pred.val[ind.1] <- 0.999
        pred.val[ind.0] <- 0.001
        rf.cv.err[j] <- sum(com != 0)/nrow(dat.test)
        rf.cv.cro[j] <- - (1/nrow(dat.test)) * sum(dat.test[,1]*log(as.numeric(pred.val)) + (1 - dat.test[,1])*log(1 - as.numeric(pred.val)))
        
      }
      
      rf.cv.errs[[i]] <- rf.cv.err
      rf.cv.cros[[i]] <- rf.cv.cro   
      
      ## EN
      
      lambda <- numeric()
      en.cv.err <- numeric()
      for (cc in 1:length(alpha)) {
        cv.kcv <- cv.glmnet(x = as.matrix(dat.trai[,-1]), y = dat.trai$y, alpha = alpha[cc], nfolds = n.folds, 
                            type.measure = "class", family = "binomial")
        lambda[cc] <- cv.kcv$lambda.min
        en.cv.err[cc] <- min(cv.kcv$cvm)
      }
      ind.min.err <- which.min(en.cv.err)[1]
      
      en.paras[[i]] <- c(alpha[ind.min.err], lambda[ind.min.err])
      en.cv.errs[[i]] <- min(en.cv.err)
      
      en.fit <- glmnet(as.matrix(dat.trai[,-1]), dat.trai$y, alpha = alpha[ind.min.err], lambda = lambda[ind.min.err], 
                       standardize = TRUE, type.measure = "class", family = "binomial")
      
      pred.val <- as.numeric(predict(en.fit, newx = as.matrix(dat.test[,-1]), type = "response"))
      pred.val <- round(pred.val, 3)
      ind.1 <- which(pred.val == 1)
      ind.0 <- which(pred.val == 0)
      pred.val[ind.1] <- 0.999
      pred.val[ind.0] <- 0.001
      (en.cv.cros[[i]] <- - (1/nrow(dat.test)) * sum(dat.test[,1]*log(as.numeric(pred.val)) + (1 - dat.test[,1])*log(1 - as.numeric(pred.val))))
      
    }
    
    nn.cv.err.dis[[q]] <- rowMeans(as.data.frame(nn.cv.errs))
    nn.cv.cro.dis[[q]] <- rowMeans(as.data.frame(nn.cv.cros))      
    
    rf.cv.err.dis[[q]] <- rowMeans(as.data.frame(rf.cv.errs))
    rf.cv.cro.dis[[q]] <- rowMeans(as.data.frame(rf.cv.cros))      
    
    en.cv.err.dis[[q]] <- rowMeans(as.data.frame(en.cv.errs))
    en.cv.cro.dis[[q]] <- rowMeans(as.data.frame(en.cv.cros))  
    
    en.para.dis[[q]] <- rowMeans(as.data.frame(en.paras))  
    
  }
  
  nn.cv.err <- as.data.frame(nn.cv.err.dis)
  nn.cv.cro <- as.data.frame(nn.cv.cro.dis)
  
  ind.opt.dis <- which.min(sapply(nn.cv.cro.dis, min))[1]
  ind.opt.para <- which.min(nn.cv.cro.dis[[ind.opt.dis]])[1]
  
  dat <- dat.list[[ind.opt.dis]]
  if (interac) {
    n.pred <- (ncol(dat)-2)/3
  }
  n.neu <- ceiling(n.pred*n.neus)
  n.opt.neu <- n.neu[ind.opt.para]
  
  nn.fit <- neuralnet(y ~ ., data = dat, hidden = c(n.opt.neu, ceiling(n.opt.neu/2), ceiling(n.opt.neu/4)), 
                      err.fct = "ce", act.fct = "logistic", algorithm = "rprop+", linear.output = FALSE)
  
  Tr.X.1 <- cbind(rep(1, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.1)[1] <- "T" 
  Tr.X.0 <- cbind(rep(0, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.0)[1] <- "T" 
  
  y.hat.1 <- neuralnet::compute(nn.fit, Tr.X.1)$net.result
  y.hat.0 <- neuralnet::compute(nn.fit, Tr.X.0)$net.result
  Z.nn <- as.numeric(y.hat.1 - y.hat.0)
  
  nn.cv.err <- nn.cv.err[ind.opt.para,]
  nn.cv.cro <- nn.cv.cro[ind.opt.para,]
  colnames(nn.cv.err) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  colnames(nn.cv.cro) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  
  out.nn <- list(cv.cro = nn.cv.cro, Z = Z.nn)
  
  rf.cv.err <- as.data.frame(rf.cv.err.dis)
  rf.cv.cro <- as.data.frame(rf.cv.cro.dis)
  
  ind.opt.dis <- which.min(sapply(rf.cv.cro.dis, min))[1]
  ind.opt.para <- which.min(rf.cv.cro.dis[[ind.opt.dis]])[1]
  
  dat <- dat.list[[ind.opt.dis]]
  if (interac) {
    n.pred <- (ncol(dat)-2)/3
  }
  mtry <- ceiling(c(n.pred*(1/2), sqrt(n.pred), log(n.pred)))
  rf.mtry <- mtry[ind.opt.para]
  
  rf.fit <- randomForest(y = as.factor(dat$y), x = dat[,-1], mtry = rf.mtry, ntree = n.trees, family = "binomial")
  
  Tr.X.1 <- cbind(rep(1, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.1)[1] <- "T" 
  Tr.X.0 <- cbind(rep(0, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.0)[1] <- "T" 
  
  y.hat.1 <- as.numeric(predict(rf.fit, Tr.X.1, "prob")[,2])
  y.hat.0 <- as.numeric(predict(rf.fit, Tr.X.0, "prob")[,2])
  Z.rf <- as.numeric(y.hat.1 - y.hat.0)
  
  rf.cv.err <- rf.cv.err[ind.opt.para,]
  rf.cv.cro <- rf.cv.cro[ind.opt.para,]
  colnames(rf.cv.err) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  colnames(rf.cv.cro) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  
  out.rf <- list(cv.cro = rf.cv.cro, Z = Z.rf, para = rf.mtry)
  
  en.cv.err <- as.data.frame(en.cv.err.dis)
  en.cv.cro <- as.data.frame(en.cv.cro.dis)
  
  ind.opt.dis <- which.min(en.cv.cro)[1]
  
  dat <- dat.list[[ind.opt.dis]]
  
  en.fit <- glmnet(as.matrix(dat[,-1]), dat$y, alpha = en.para.dis[[ind.opt.dis]][1], lambda = en.para.dis[[ind.opt.dis]][2], 
                   standardize = TRUE, type.measure = "class", family = "binomial")
  
  Tr.X.1 <- cbind(rep(1, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.1)[1] <- "T" 
  Tr.X.0 <- cbind(rep(0, nrow(dat)), dat[,-c(1,2)])
  colnames(Tr.X.0)[1] <- "T" 
  
  y.hat.1 <- as.numeric(predict(en.fit, newx = as.matrix(Tr.X.1), type = "response"))
  y.hat.0 <- as.numeric(predict(en.fit, newx = as.matrix(Tr.X.0), type = "response"))
  Z.en <- as.numeric(y.hat.1 - y.hat.0)
  
  colnames(en.cv.err) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  colnames(en.cv.cro) <- c("Euclidean", "Jaccard", "BC", "UUniFrac", "GUniFrac", "WUniFrac")
  
  out.en <- list(cv.cro = en.cv.cro, Z = Z.en, para = en.para.dis[[ind.opt.dis]])
  
  ind.opt.1 <- which.min(c(min(nn.cv.cro), min(rf.cv.cro), min(en.cv.cro)))[1]
  
  if (ind.opt.1 == 1) {
    out.dml <- out.nn$Z
  }
  
  if (ind.opt.1 == 2) {
    out.dml <- out.rf$Z
  }
  
  if (ind.opt.1 == 3) {
    out.dml <- out.en$Z
  }
  
  return(list(out.en = out.en, out.rf = out.rf, out.dfn = out.nn, Z = out.dml))
  
}
