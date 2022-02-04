\name{biom.qc}
\alias{biom.qc}
\title{
Quality Controls and Data Transformation
}
\description{
This function performs quality controls and data transformations that are needed for MiVT.
}
\usage{
biom.qc(biom, kingdom = "Bacteria", lib.size.cut.off = 1000, mean.prop.cut.off = 0, rem.tax.com = c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA"), rem.tax.par = c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown"))
}

\arguments{
  \item{biom}{
A microbiome data in the phyloseq format. sample_data(biom) should contain two binary variables: y (response) and Tr (treatment). See the example data using gen.syn.dat().
}
  \item{kingdom}{
A microbial kingdom to be analyzed, such as 'Bacteria', 'Archaea', 'Eukaryota' or 'all'. 'all' is for all kingdoms in the taxonomic table (Default: 'Bacteria').
}
  \item{lib.size.cut.off}{
A minimum total read count for subjects to keep in the microbiome data (Default: 1000).
}
  \item{mean.prop.cut.off}{
A minimum mean proportion for microbial features (OTUs or ASVs) to keep in the microbiome data (Default: 0).
}
  \item{rem.tax.com}{
Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings (Default: c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA")).
}
  \item{rem.tax.par}{
Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (Default: c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown")).
}
}

\value{
$tax.prop: A list of tables for the proportions of microbial taxa on each taxonomic rank (Phylum, Class, Order, Family, Genus, Species).
$otu.tab: A feature (OTU or ASV) table where rows are features and columns are subjects.
$tax.tab: A taxonomic table where rows are features (OTUs or ASVs), and columns are seven taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species).
$sam.dat: A metadata/sample information where rows are subjects and columns are variables. It should contain two binary variables: y (response) and Tr (treatment).
$tree: A rooted phylogenetic tree.
}

\references{
Koh, H. Subgroup identification using virtual twins for human microbiome studies. (Under review).
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
}
