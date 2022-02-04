biom.clean <-
function(biom, kingdom = "all", lib.size.cut.off = 1000, mean.prop.cut.off = 0,
                       rem.tax.com = c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA"),
                       rem.tax.par = c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown")) {
  
  tax.tab <- tax_table(biom)
  if (kingdom != "all") {
    ind <- is.element(tax.tab[,1], kingdom)
    rownames(tax.tab)[ind]
    biom <- prune_taxa(rownames(tax.tab)[ind], biom)
  }
  
  otu.tab <- otu_table(biom)
  tax.tab <- tax_table(biom)
  tree <- phy_tree(biom)
  sam.dat <- sample_data(biom)
  
  ind.otu.sam <- is.element(dim(otu.tab), nrow(sam.dat))
  if (sum(ind.otu.sam) == 0) {
    stop("the numbers of subjects (n) in OTU table and sample data are not the same")
  }
  if (sum(ind.otu.sam) == 1 & ind.otu.sam[1]) {
    otu.tab <- t(otu.tab)
  }
  if (!identical(colnames(otu.tab), rownames(sam.dat))) {
    stop("subject IDs (colunm names of OTU table and row names of sample data) are not the same")
  }
  ind.com.otu <- intersect(intersect(rownames(otu.tab), tree$tip.label), rownames(tax.tab))
  if (length(ind.com.otu) == 0) {
    stop("there is no common OTUs among OTU table, taxonomic table and tree tip labels")
  }
  
  ind.com.1 <- which(rownames(otu.tab) %in% ind.com.otu)
  otu.tab <- otu.tab[ind.com.1,]
  ind.com.2 <- which(rownames(tax.tab) %in% ind.com.otu)
  tax.tab <- tax.tab[ind.com.2,]
  tree <- prune_taxa(ind.com.otu, tree)
  if(!is.rooted(tree)) {
    tree <- phangorn::midpoint(tree)
  }
  tax.tab <- tax.tab.clean(tax.tab, rem.tax.com = rem.tax.com, rem.tax.par = rem.tax.par)
  biom <- merge_phyloseq(otu.tab, tax.tab, tree, sam.dat)
  biom <- otu.tab.clean(biom, lib.size.cut.off, mean.prop.cut.off)
  
  return(biom)  
}
