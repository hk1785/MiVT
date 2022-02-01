boot.test <-
function(Z.hat, tree.fit, n.boot = 20000) {
  
  if (length(unique(tree.fit$where)) >= 1) {
    out <- as.data.frame(c(length(tree.fit$where), mean(Z.hat), mean(Z.hat), 0, 1))
    rownames(out) <- c("N", "Overall TE", "Subgroup TE", "Subgroup TE - Overall TE", "P-value (BoRT)")
    colnames(out) <- "Overall(No partition)"
  }  
  
  if (length(unique(tree.fit$where)) >= 2) {
    leaves <- names(table(tree.fit$where))
    ns <- numeric()
    ove.trt.effs <- numeric()
    sub.trt.effs <- numeric()
    Q.ests <- numeric()
    pvs <- numeric()
    for (k in 1:length(leaves)) {
      ind.sub <- which(tree.fit$where == leaves[k])
      ns[k] <- n <- length(ind.sub)
      ove.trt.effs[k] <- ove.trt.eff <- mean(Z.hat)
      sub.trt.effs[k] <- sub.trt.eff <- mean(Z.hat[ind.sub])
      Q.ests[k] <- Q.est <- sub.trt.eff - ove.trt.eff
      Q.est.null <- numeric()
      for (j in 1:n.boot) {
        Z.hat.res <- sample(Z.hat, replace = TRUE)
        Q.est.null[j] <- mean(Z.hat.res[ind.sub]) - mean(Z.hat.res)
      }
      pvs[k] <- (sum(abs(Q.est.null) >= abs(Q.est)) + 1)/(n.boot + 1)
    }
    out <- as.data.frame(rbind(ns, ove.trt.effs, sub.trt.effs, Q.ests, pvs))
    rownames(out) <- c("N", "Overall TE", "Subgroup TE", "Subgroup TE - Overall TE", "P-value (BoRT)")
    colnames(out) <- leaves
  }
  
  return(out)
}
