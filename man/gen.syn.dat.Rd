\name{gen.syn.dat}
\alias{gen.syn.dat}
\title{
Generate Example Microbiome Data
}
\description{
This function generates example synthetic data that are required to use MiVT.
}
\usage{
gen.syn.dat(tree, tax.tab, prop, disp, num.sams = 50, seq.depth = sample(10000:1e+05, 50), keep.cut.off = 200)
}
\arguments{
  \item{tree}{
A rooted phylogenetic tree.
}
  \item{tax.tab}{
A taxonomic table where rows are features (OTUs or ASVs), and columns are seven taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species)
}
  \item{prop}{
A vector of proportion parameters for the Dirichlet-Multinomial distribution.
}
  \item{disp}{
A dispersion parameter for the Dirichlet-Multinomial distribution.
}
  \item{num.sams}{
A sample size (Default: 50).
}
  \item{seq.depth}{
A vector of total read counts across subjects (Default: sample(10000:1e+05, 50)).
}
  \item{keep.cut.off}{
A number of features to keep in the microbiome data (=< 755) (Default: 200).
}
}

\value{
A synthetic microbiome data in phyloseq format.
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
}
