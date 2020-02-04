# receptLoss

This package contains functions to detect nuclear hormon receptor (NHR) expression loss in subsets of tumors relative to normal tissue:

* `receptLoss`: the main function in this package.
Input: 2 matrices of expression data, one from tumor and one from adj normal tissue.
Output: Table of summary statistics, with each row as a gene and each column a separate statistic. The delta mu value is

* `plotReceptLoss`: A function to visualize the output from receptLoss for a particular gene.

It also includes a dataset consisting of all NHRs:

* `nhrs` (Source: http://www.guidetopharmacology.org/DATA/targets_and_families.csv)

## Installation

```R
# Install the development version from GitHub
devtools::install_github("dpique/receptLoss")
```


