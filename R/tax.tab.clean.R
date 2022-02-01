tax.tab.clean <-
function(tax.tab, rem.tax.com = rem.tax.com, rem.tax.par = rem.tax.par, na.code = "NANANA") {
  
  tax.tab.c <- tax.tab
  for (i in 1:ncol(tax.tab)) {
    taxa <- as.character(tax.tab.c[,i])
    tax.tab.c[is.na(taxa), i] <- na.code
    tax.tab.c[is.element(taxa, rem.tax.com), i] <- na.code
    tax.tab.c[grep(paste(rem.tax.par, collapse="|"), taxa), i] <- na.code
    uniq.taxa <- names(table(taxa))
    for (j in 1:length(uniq.taxa)) {
      tax.tab.c[is.element(taxa, paste(uniq.taxa[j], 1:100)), i] <- uniq.taxa[j]
    }
  }
  
  return(tax.tab.c)
}
