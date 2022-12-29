# MiVT

Title: Microbiome Virtual Twins

Version: 1.0

Date: 2022-12-05

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hyunwook.koh@stonybrook.edu>

Description: This R package provides facilities for MiVT that predicts treatment effects using dML, and then identifies subgroups by treatment effects based on microbiome composition using BoRT to evaluate the interplay between microbiome and treatment. 

NeedsCompilation: no

Depends: R(>= 4.1.1), cluster, compositions, dirmult, glmnet, GUniFrac, ecodist, neuralnet, phangorn, phyloseq, proxy, randomForest, rpart, rpart.plot, splitTools, zCompositions

License: GPL-2

NeedsCompilation: no

URL: https://github.com/hk1785/MiVT

## Reference

* Koh, H. Subgroup identification using virtual twins for human microbiome studies. *_Under review_*

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/MiVT/issues) or email Hyunwook Koh (hyunwook.koh@stonybrook.edu).

## Prerequites

Rtools
```
https://cran.r-project.org/bin/windows/Rtools/rtools40.html
```
phyloseq
```
https://joey711.github.io/phyloseq/
```
cluster, compositions, devtools, dirmult, glmnet, GUniFrac, ecodist, neuralnet, phangorn, proxy, randomForest, rpart, rpart.plot, splitTools, zCompositions
```
install.packages(c("cluster", "compositions", "devtools", "dirmult", "glmnet", "GUniFrac", "ecodist", "nuralnet", "phangorn", "proxy", "randomForest", "rpart", "rpart.plot", "splitTools", "zCompositions"))
```

## Installation

```
library(devtools)
install_github("hk1785/MiVT", force=T)
```

---------------------------------------------------------------------------------------------------------------------------------------

# Manual
This R package <MiVT> contains four core functions, gen.syn.dat, biom.qc, dML and BoRT. Details are below.

## :mag: gen.syn.dat

### Description
This function generates example microbiome data.

### Usage
```
gen.syn.dat(tree, tax.tab, prop, disp, num.sams = 50, seq.depth = sample(10000:1e+05, 50), keep.cut.off = 200)
```

### Arguments
* _tree_ - A rooted phylogenetic tree.
* _tax.tab_ - A taxonomic table where rows are features (OTUs or ASVs), and columns are seven taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species)
* _prop_ - A vector of proportion parameters for the Dirichlet-Multinomial distribution.\
* _disp_ - A dispersion parameter for the Dirichlet-Multinomial distribution.
* _num.sams_ - A sample size (Default: 50).
* _seq.depth_ - A vector of total read counts across subjects (Default: sample(10000:1e+05, 50)).
* _keep.cut.off_ - A number of features to keep in the microbiome data (=< 755) (Default: 200).

### Values
A synthetic microbiome data in the 'phyloseq' format.

### Example
Import requisite R packages
```
library('phyloseq')
library('cluster')
library('dirmult')
library('phangorn')
library('compositions')
library('zCompositions')
library('GUniFrac')
library('ecodist') 
library('proxy')
library('glmnet')
library('randomForest')
library('neuralnet')
library('splitTools')
library('rpart')
library('rpart.plot')
library('MiVT')
```
Generate example microbiome data
```
data(fit)
data(tree)
data(tax.tab)

prop <- fit$pi
disp <- fit$theta

sim.biom <- gen.syn.dat(tree = tree, tax.tab = tax.tab, prop = prop, disp = disp)
sim.biom
```

## :mag: biom.qc

### Description
This function performs quality controls and data transformations that are needed for MiVT.

### Usage
```
biom.qc(biom = biom, kingdom = "Bacteria", lib.size.cut.off = 1000, mean.prop.cut.off = 0, rem.tax.com = c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA"), rem.tax.par = c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown"))
```

### Arguments
* _biom_ - A microbiome data in the phyloseq format. sample_data(biom) should contain two binary variables: y (response) and Tr (treatment). See the example data using gen.syn.dat().
* _kingdom_ - A microbial kingdom to be analyzed, such as 'Bacteria', 'Archaea', 'Eukaryota' or 'all'. 'all' is for all kingdoms in the taxonomic table (Default: 'Bacteria').
* _lib.size.cut.off_ - A minimum total read count for subjects to keep in the microbiome data (Default: 1000).
* _mean.prop.cut.off_ - A minimum mean proportion for microbial features (OTUs or ASVs) to keep in the microbiome data (Default: 0).
* _rem.tax.com_ - Remove taxonomic names in the taxonomic table that are completely matched with the specified character strings (Default: c("", "gut metagenome", "mouse gut metagenome", "metagenome", "NANANA")).
* _rem.tax.par_ - Remove taxonomic names in the taxonomic table that are partially matched with the specified character strings (Default: c("uncultured", "incertae", "Incertae", "unclassified", "unidentified", "unknown")).

### Values
$tax.prop: A list of tables for the proportions of microbial taxa on each taxonomic rank (Phylum, Class, Order, Family, Genus, Species).

$otu.tab: A feature (OTU or ASV) table where rows are features and columns are subjects.

$tax.tab: A taxonomic table where rows are features (OTUs or ASVs), and columns are seven taxonomic ranks (Kingdom, Phylum, Class, Order, Family, Genus, Species).

$sam.dat: A metadata/sample information where rows are subjects and columns are variables. It should contain two binary variables: y (response) and Tr (treatment).

$tree: A rooted phylogenetic tree.

### References
* Koh, H. Subgroup identification using virtual twins for human microbiome studies. *_Under review_*

### Example
Import requisite R packages
```
library('phyloseq')
library('cluster')
library('dirmult')
library('phangorn')
library('compositions')
library('zCompositions')
library('GUniFrac')
library('ecodist') 
library('proxy')
library('glmnet')
library('randomForest')
library('neuralnet')
library('splitTools')
library('rpart')
library('rpart.plot')
library('MiVT')
```
Generate example microbiome data
```
data(fit)
data(tree)
data(tax.tab)

prop <- fit$pi
disp <- fit$theta

sim.biom <- gen.syn.dat(tree = tree, tax.tab = tax.tab, prop = prop, disp = disp)
sim.biom
```
Perform quality controls and data transformations
```
qc.out <- biom.qc(biom = sim.biom)
```

## :mag: dML

### Description
This function implements dML to predicts treatment effects.

### Usage
```
dML(y, Tr, X, tree, n.folds = 10, n.rep = 2, alpha = seq(0.05, 0.95, 0.05), n.trees = 1000, n.neus = c(1/2, 1/3, 1/4))
```

### Arguments
* _y_ - A vector of binary responses.
* _Tr_ - A vector of binary treatment status. 
* _X_ - A feature (OTU or ASV) table where rows are features and columns are subjects.
* _tree_ - A rooted phylogenetic tree.
* _n.folds_ - The number of folds in the k-fold cross-validation (Default: 10).
* _n.rep_ - The number of repeats of the k-fold cross-validation (Default: 2).
* _alpha_ - A rooted phylogenetic tree.
* _n.trees_ - The number of folds in the k-fold cross-validation (Default: 10).
* _n.neus_ - Candidate numbers of neurons on the first hidden layer for the three layer deep feedforward network: {[M/2], [M/4], [M/8]}, {[M/3], [M/6], [M/12]}, and {[M/4], [M/8], [M/16]}, where M is the number of predictors (coordinates). 

### Values
$out.en$cv.cro: CV cross-entropy values for the elastic net and each distance measure.
$out.en$Z: Predicted treatment effects using the elastic net.

$out.rf$cv.cro: CV cross-entropy values for the random forest and each distance measure.
$out.rf$Z: Predicted treatment effects using the random forest.

$out.dfn$cv.cro: CV cross-entropy values for the deep feedforward network and each distance measure.
$out.dfn$Z: Predicted treatment effects using the deep feedforward network.

$Z: Predicted treatment effects using dML.

### References
* Koh, H. Subgroup identification using virtual twins for human microbiome studies. *_Under review_*

* Foster, J. C., Taylor, J. M. & Ruberg, S. J. Subgroup identification from randomized clinical trial data. *_Stat. Med._* 30(24), 2867-2880 (2011).

### Example
Import requisite R packages
```
library('phyloseq')
library('cluster')
library('dirmult')
library('phangorn')
library('compositions')
library('zCompositions')
library('GUniFrac')
library('ecodist') 
library('proxy')
library('glmnet')
library('randomForest')
library('neuralnet')
library('splitTools')
library('rpart')
library('rpart.plot')
library('MiVT')
```
Generate example microbiome data
```
data(fit)
data(tree)
data(tax.tab)

prop <- fit$pi
disp <- fit$theta

sim.biom <- gen.syn.dat(tree = tree, tax.tab = tax.tab, prop = prop, disp = disp)
sim.biom
```
Perform quality controls and data transformations
```
qc.out <- biom.qc(biom = sim.biom)
```
Perform dML
```
dml.out <- dML(y = qc.out$sam.dat$y, Tr = qc.out$sam.dat$Tr, X = qc.out$otu.tab, tree = qc.out$tree)
```

## :mag: BoRT

### Description
This function implements BoRT for subgroup identification and significance testing.

### Usage
```
BoRT(Z, tax.prop, tax.rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species"), minsplit = 10, minbucket = 5, cp = 0.01, n.boot = 20000)
```

### Arguments
* _Z_ - Predicted treatment effects using dML (see dML()).
* _tax.prop_ - A list of tables for the proportions of microbial taxa on each taxonomic rank (Phylum, Class, Order, Family, Genus, Species) (see biom.clean()).
* _tax.rank_ -A taxonomic rank to be used for the subgroup identification (tax.rank = c("Phylum", "Class", "Order", "Family", "Genus", "Species")).
* _minsplit_ - The minimum number of observations that must exist in a node in order for a split to be attempted (Default = 10).
* _minbucket_ - The minimum number of observations in any terminal node (Default = 5).
* _cp_ - The complexity parameter of the decision tree (Default = 0.01).
* _n.boot_ - The number of bootstrap samples (Default = 20000).

### Values
$Sel.Taxa: Short taxonomic IDs and full taxonomic names.

$BoRT.out: The output table of BoRT. Columns are the identified subgroups that correspond with the terminal nodes from left to right. N is the sample size for each subgroup. Overall TE represents the overall treatment effect, and Subgroup TE represents the subgroup treatment effect. 


### References
* Koh, H. Subgroup identification using virtual twins for human microbiome studies. *_Under review_*

* Foster, J. C., Taylor, J. M. & Ruberg, S. J. Subgroup identification from randomized clinical trial data. *_Stat. Med._* 30(24), 2867-2880 (2011).

### Example
Import requisite R packages
```
library('phyloseq')
library('cluster')
library('dirmult')
library('phangorn')
library('compositions')
library('zCompositions')
library('GUniFrac')
library('ecodist') 
library('proxy')
library('glmnet')
library('randomForest')
library('neuralnet')
library('splitTools')
library('rpart')
library('rpart.plot')
library('MiVT')
```
Generate example microbiome data
```
data(fit)
data(tree)
data(tax.tab)

prop <- fit$pi
disp <- fit$theta

sim.biom <- gen.syn.dat(tree = tree, tax.tab = tax.tab, prop = prop, disp = disp)
sim.biom
```
Perform quality controls and data transformations
```
qc.out <- biom.qc(biom = sim.biom)
```
Perform dML
```
dml.out <- dML(y = qc.out$sam.dat$y, Tr = qc.out$sam.dat$Tr, X = qc.out$otu.tab, tree = qc.out$tree)
```
Perform BoRT
```
bort.out <- BoRT(Z = dml.out$Z, tax.prop = qc.out$tax.prop, tax.rank = "Genus")
bort.out
```
