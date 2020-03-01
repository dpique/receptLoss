receptLoss v0.99.8 (2020-03-01)
==============

* Updated bug in 'toMatrix' function in scripts.R, allowing it to accept the
objects with multiple classes


receptLoss v0.99.6 (2020-03-01)
==============

* Updated the 'toMatrix' function in scripts.R, allowing it to accept the
RangedSummarizedExperiment class, and added option for specifying rownames

receptLoss v0.99.5 (2020-02-14)
==============

* update .Rbuildignore


receptLoss v0.99.4 (2020-02-12)
==============

* removed .Rproj files from git repo

* Updated gitignore file to include .Rproj files

receptLoss v0.99.3 (2020-02-12)
==============

* removed .DS_Store and .Rproj files from git repo


receptLoss v0.99.2 (2020-02-12)
==============

* exported 'toMatrix' function to fix R CMD check error.


receptLoss v0.99.1 (2020-02-11)
==============

* As per Martin Morgan's feedback from BioConductor review, 
added a function, 'toMatrix', to accept SummarizedExperiment 
and other matrix-like objects


receptLoss v0.99.0 (2020-02-05)
==============

Changes:

* Passed R CMD Check and devtools::build()
* Submitted to Bioconductor
