\name{GGPA2-package}
\alias{GGPA2-package}
\docType{package}
\title{
graph-GPA 2.0
}
\description{
This package provides functions for fitting graph-GPA 2.0, a graphical model for multi-disease analysis of GWAS results with integration of functional annotation data
}
\details{
\tabular{ll}{
Package: \tab GGPA2\cr
Type: \tab Package\cr
Version: \tab 1.9.1\cr
Date: \tab 2022-04-10\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains a main class, \code{GGPA2}, which represents graph-GPA 2.0 model fit.
This package contains four main methods, 
\code{GGPA2}, \code{assoc}, and \code{plot}.
\code{GGPA2} method fits the graph-GPA 2.0 model
and \code{assoc} method implements association mapping.
\code{plot} method provides a graph representing genetic relationship among phenotypes.
}
\author{
Hang J. Kim, Dongjun Chung, and Qiaolan Deng

Maintainer: Hang J. Kim <hang.kim@uc.edu>, Dongjun Chung <chung.911@osu.edu>
}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388
}
\keyword{ package }
\seealso{
\code{\link{GGPA2}}, \code{\link{assoc}}, \code{\link{plot}}, 
\code{\linkS4class{GGPA2}}.
}
\examples{
# load simulation data

data(simulationAnnot)

# fit graph-GPA 2.0 model

fit <- GGPA2( simulationAnnot$pmat, nBurnin=200, nMain=200 )
fit

# fit graph-GPA 2.0 model using a prior phenotype graph
# as an example, edge 6-7 added & edge 2-3 removed in pgraph

pgraph <- matrix( 0, ncol(simulationAnnot$pmat), ncol(simulationAnnot$pmat) )
pgraph[1,2] <- pgraph[1,3] <- pgraph[6,7] <- pgraph[4,5] <- 1
fit_prior <- GGPA2( simulationAnnot$pmat, pgraph = pgraph, nBurnin=200, nMain=200 )
fit_prior


# fit graph-GPA 2.0 model using a prior phenotype graph and functional annotations

fit_annot <- GGPA2( simulationAnnot$pmat, pgraph = pgraph, annotMat = simulationAnnot$A_mat, nBurnin=200, nMain=200 )


# association mapping for each phenotype

head(assoc( fit_annot, FDR=0.1, fdrControl="global" ))

# hypothesis testing for 1st and 2nd phenotype pair

head(assoc( fit_annot, FDR=0.1, fdrControl="global", i=1, j=2 ))

# plot phenotype graph

plot(fit)
plot(fit_prior)
plot(fit_annot)




}
