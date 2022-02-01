gen.syn.dat <-
function(tree, tax.tab, prop, disp, 
                        num.sams = 50, seq.depth = sample(10000:100000, 50), keep.cut.off = 200) {
  
  prop.sig.X = 0.05; a0 = 0.1; a2 = 0; b0 = 0.5; b2 = 0.1; b1 = 0.05

  sim.otu.tab <- simPop(J = num.sams, n = seq.depth, pi = prop, theta = disp)$data*(100/seq.depth)
  colnames(sim.otu.tab) <- paste("X", 1:length(tree$tip.label), sep = "")
  
  keep.otus <- names(colSums(sim.otu.tab)[order(colSums(sim.otu.tab), decreasing = TRUE)])[1:keep.cut.off]
  ind <- match(keep.otus, colnames(sim.otu.tab))
  X.ori <- X <- sim.otu.tab[,ind]
  sim.tree <- prune_taxa(colnames(X), tree)
  
  ind <- match(sim.tree$tip.label, colnames(X))
  X.ori <- X <- X[,ind]
  
  ind <- match(keep.otus, rownames(tax.tab))
  sim.tax.tab <- tax.tab[ind,]
  
  ind <- match(sim.tree$tip.label, rownames(sim.tax.tab))
  sim.tax.tab <- sim.tax.tab[ind,]
  
  identical(colnames(X), rownames(sim.tax.tab))
  
  Tr <- c(rep(0, num.sams/2), rep(1, num.sams/2))
  
  T0X <- abs(Tr-1)*X
  colnames(T0X) <- paste("T0", colnames(X), sep = "")
  
  T1X <- Tr*X
  colnames(T1X) <- paste("T1", colnames(X), sep = "")
  
  Tr.X <- cbind(Tr, X, T0X, T1X)
  colnames(Tr.X)[1] <- "T" 
  
  # Predictors
  
  ind.pl <- which(Tr == 0)    
  ind.tr <- which(Tr == 1)
  
  X.pl <- X[ind.pl,]
  X.tr <- X[ind.tr,]
  
  pam.out <- pam(cophenetic(sim.tree), round(1/prop.sig.X))
  mean.abu <- numeric()
  for (i in 1:round(1/prop.sig.X)) {
    sel.taxa <- names(which(pam.out$clustering == i))
    mean.abu[i] <- mean(rowSums(X[,which(colnames(X) %in% sel.taxa)]))
  }
  ind.can <- intersect(which(mean.abu >= 5), which(mean.abu < 15))
  sel.clust <- sample(ind.can, 1)
  sel.taxa <- names(which(pam.out$clustering == sel.clust))
  sig.X.tr.ind  <- sig.X.pl.ind <- which(colnames(X) %in% sel.taxa)
  
  # Responses
  
  beta1.pl <- rep(a0, length(ind.pl))  
  beta2.pl <- rep(a2, length(sig.X.pl.ind))
  beta1.tr <- rep(b0, length(ind.tr))
  beta2.tr <- rep(b2, length(sig.X.tr.ind))  
  beta2 <- rep(b1, length(sig.X.tr.ind))
  
  ind.pl <- which(Tr == 0)    
  ind.tr <- which(Tr == 1)
  
  res.pl <- beta1.pl + X.pl[,sig.X.tr.ind] %*% beta2 + X.pl[,sig.X.pl.ind] %*% beta2.pl
  res.p.pl <- 1/(1 + exp(-res.pl))
  y.pl <- rbinom(length(res.p.pl), 1, res.p.pl)
  
  res.tr <- beta1.tr + X.tr[,sig.X.tr.ind] %*% beta2 + X.tr[,sig.X.tr.ind] %*% beta2.tr
  res.p.tr <- 1/(1 + exp(-res.tr))
  y.tr <- rbinom(length(res.p.tr), 1, res.p.tr)
  
  y <- c(y.pl, y.tr)
  X <- X*(seq.depth/100)
  
  sim.sam.dat <- as.data.frame(cbind(y, Tr))
  rownames(sim.sam.dat) <- paste("P", 1:num.sams, sep = "")
  sim.sam.dat <- sample_data(sim.sam.dat)
  rownames(X) <- paste("P", 1:num.sams, sep = "")
  sim.otu.tab <- otu_table(t(X), taxa_are_rows = TRUE)
  sim.tax.tab <- tax_table(sim.tax.tab)
  
  identical(rownames(sim.tax.tab), rownames(sim.otu.tab))
  identical(rownames(sim.tax.tab), sim.tree$tip.label)
  identical(rownames(sim.otu.tab), sim.tree$tip.label)
  identical(colnames(sim.otu.tab), rownames(sim.sam.dat))

  sim.biom <- merge_phyloseq(sim.otu.tab, sim.tax.tab, sim.sam.dat, sim.tree)
  
  return(sim.biom)
  
}
