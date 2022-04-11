graph-GPA 2.0
===

graph-GPA 2.0 is a graphical model for prioritizing GWAS results and investigating pleiotropic architecture and enrichment for functional annotation within a unified framework (Deng et al., 2022; Chung et al., 2017). 'GGPA2' package provides user-friendly interface to fit graph-GPA 2.0 models, implement association mapping, generate a phenotype graph, and implement enrichment analysis for functional annotation. It also allows to utilize a prior phenotype graph (Kim et al., 2017).

DDNet
===

DDNet is a web interface to infer the disease-disease network from the literature mining. Using DDNet, users can query diseases of interest, investigate relationship among these diseases visually and dynamically, and download various information including an adjacency matrix. This adjacency matrix can be used as a prior phenotype graph to guide the estimation of pleiotropic architecture. Please refer Kim et al. (2017) and the R package vignette for more details. DDNet is accessible from the following web address: 

http://chunglab.io/ddnet/

Installation
===========

To install the development version of GGPA2, it's easiest to use the 'devtools' package. Note that GGPA depends on the 'Rcpp' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.

```
#install.packages("devtools")
library(devtools)
install_github("dongjunchung/GGPA2")
```

Usage
===========

The R package vignette will provide a good start point for the genetic analysis using GGPA2 package, including the overview of GGPA2 package and the example command lines:

```
library(GGPA2)
vignette("GGPA2-example")
```
The following two help pages will also provide quick references for GGPA2 package and the example command lines:

```
package?GGPA2
class?GGPA2
```

References
==========

Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139–2141.

Chung D, Kim H, and Zhao H (2017), "graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," *PLOS Computational Biology*, 13(2): e1005388.
