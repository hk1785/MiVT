biom.qc <-
function(biom = biom, kingdom = "Bacteria", lib.size.cut.off = 1000, mean.prop.cut.off = 0, 
                    rem.tax.com = c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA"),
                    rem.tax.par = c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown")) {
  
  biom.after.qc <- biom.clean(biom = biom, kingdom = kingdom, lib.size.cut.off = lib.size.cut.off, mean.prop.cut.off = mean.prop.cut.off, rem.tax.com = rem.tax.com, rem.tax.par = rem.tax.par)
  
  qc.otu.tab <- otu_table(biom.after.qc)
  qc.tax.tab <- tax_table(biom.after.qc)
  qc.sam.dat <- sample_data(biom.after.qc)
  qc.tree <- phy_tree(biom.after.qc)
  
  rare.biom <- rarefy.func(biom.after.qc, cut.off = min(colSums(qc.otu.tab)), multi.rarefy = 1)
  
  ra.otu.tab <- otu_table(rare.biom)
  ra.tax.tab <- tax_table(rare.biom)
  ra.sam.dat <- sample_data(rare.biom)
  ra.tree <- phy_tree(rare.biom)
  
  taxa.out <- tax.trans(qc.otu.tab, qc.tax.tab, ra.otu.tab, ra.tax.tab)
  tax.prop <- taxa.out$prop
  
  out <- list(tax.prop = tax.prop, otu.tab = qc.otu.tab, tax.tab = qc.tax.tab, sam.dat = qc.sam.dat, tree = qc.tree)
  
  return(out)
}
