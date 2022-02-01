\name{BoRT}
\alias{BoRT}
\title{
Bootstrap-based Test Regression Tree
}
\description{
This function implements BoRT. 
}
\usage{
BoRT(Z, tax.prop, tax.rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), minsplit = 10, minbucket = 5, cp = 0.01, n.boot = 20000)
}
\arguments{
  \item{Z}{
Predicted treatment effects using dML (see dML()).
}
  \item{tax.prop}{
A list of tables for the proportions of microbial taxa on each taxonomic rank (Phylum, Class, Order, Family, Genus, Species) (see biom.clean()).
}
  \item{tax.rank}{
A taxonomic rank to be used for the subgroup identification (tax.rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species")).
}
  \item{minsplit}{
The minimum number of observations that must exist in a node in order for a split to be attempted (Default = 10).
}
  \item{minbucket}{
The minimum number of observations in any terminal node (Default = 5).
}
  \item{cp}{
The complexity parameter of the decision tree (Default = 0.01).
}
  \item{n.boot}{
The number of bootstrap samples (Default = 20000).
}
}

\value{
$Sel.Taxa: Short taxonomic IDs and full taxonomic names.
$BoRT.out: The output table of BoRT. Columns are the identified subgroups that correspond with the terminal nodes from left to right. N is the sample size for each subgroup. Overall TE represents the overall treatment effect, and Subgroup TE represents the subgroup treatment effect. 
}

\references{
Koh, H. Subgroup identification using virtual twins for human microbiome studies. (Under review).

Foster, J. C., Taylor, J. M. & Ruberg, S. J. Subgroup identification from randomized clinical trial data. Stat. Med. 30(24), 2867-2880 (2011).
}

\author{
Hyunwook Koh
}

\examples{
data(fit)
data(tree)
data(tax.tab)

prop <- fit$pi
disp <- fit$theta

sim.biom <- gen.syn.dat(tree = tree, tax.tab = tax.tab, prop = prop, disp = disp)
sim.biom

qc.out <- biom.qc(biom = sim.biom)

dml.out <- dML(y = qc.out$sam.dat$y, Tr = qc.out$sam.dat$Tr, X = qc.out$otu.tab, tree = qc.out$tree)

bort.out <- BoRT(Z = dml.out$Z, tax.prop = qc.out$tax.prop, tax.rank = "Genus")
bort.out


}