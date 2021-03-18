# MOVICS 0.99.5

* Minimal version of R was declined to 4.0.1.
* Two functions were added. 
* Troubleshooting section was updated.

# MOVICS 0.99.6

* Add an error notification prompt: current version of MOVICS can support up to 6 omics datasets and at least 2 datasets must be provided.
* Fix bugs when providing more than one binary matrix.
* Remove argument of `is.binary` in `getIntNMF()`, please use `type` instead.

# MOVICS 0.99.7

* Fix bugs when setting argument of `N.clust` as 2 in `getConsensusClustering()`.

# MOVICS 0.99.8

* Fix bugs where argument of `fig.path` is useless in `getSilhouette()`.

# MOVICS 0.99.9

* Fix bugs in `compTMB()` if the maximum number of x-axis is less than 3.

# MOVICS 0.99.10

* Fix bugs in `compAgree()` if argument of `subt2comp` has more than 2 variables.
* Move heatmap legend and annotation legend to the left in `runMarker()` and `runGSVA()`, and return the complexheatmap object for further modification if necessary.

# MOVICS 0.99.11

* Fix bugs in `compFGA()` that poses error when only 2 subtypes were identified.
