\name{dML}
\alias{dML}
\title{
Distance-based Machine Learning
}
\description{
This function implements dML to predicts treatment effects.
}
\usage{
dML(y, Tr, X, tree, n.folds = 10, n.rep = 2, alpha = seq(0.05, 0.95, 0.05), n.trees = 1000, n.neus = c(1/2, 1/3, 1/4))
}

\arguments{
  \item{y}{
A vector of binary responses.
}
  \item{Tr}{
A vector of binary treatment status. 
}
  \item{X}{
A feature (OTU or ASV) table where rows are features and columns are subjects.
}
  \item{tree}{
A rooted phylogenetic tree.
}
  \item{n.folds}{
The number of folds in the k-fold cross-validation (Default: 10).
}
  \item{n.rep}{
The number of repeats of the k-fold cross-validation (Default: 2).
}
  \item{alpha}{
Candidate alpha values that modulates the amount of L1 and L2 penalties to their linear combination for the elastic net (Default: seq(0.05, 0.95, 0.05)).
}
  \item{n.trees}{
The number of bagged decision trees for the random forest (Default: 1000).
}
  \item{n.neus}{
Candidate numbers of neurons on the first hidden layer for the three layer deep feedforward network: {[M/2], [M/4], [M/8]}, {[M/3], [M/6], [M/12]}, and {[M/4], [M/8], [M/16]}, where M is the number of predictors (coordinates). 
}
}

\value{
$out.en$cv.cro: CV cross-entropy values for the elastic net and each distance measure.
$out.en$Z: Predicted treatment effects using the elastic net.

$out.rf$cv.cro: CV cross-entropy values for the random forest and each distance measure.
$out.rf$Z: Predicted treatment effects using the random forest.

$out.dfn$cv.cro: CV cross-entropy values for the deep feedforward network and each distance measure.
$out.dfn$Z: Predicted treatment effects using the deep feedforward network.

$Z: Predicted treatment effects using dML.
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
}