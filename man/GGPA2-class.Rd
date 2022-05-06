\name{GGPA2-class}
\alias{GGPA2-class}
\alias{show,GGPA2-method}
\alias{plot,GGPA2-method}
\alias{plotAnnot,GGPA2-method}
\alias{plot,GGPA2,missing-method}
\alias{fdr,GGPA2-method}
\alias{estimates,GGPA2-method}
\alias{fdr}
\alias{estimates}
\Rdversion{1.1}
\docType{class}
\title{Class "GGPA2" }
\description{
This class represents graph-GPA 2.0 model fit.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GGPA2", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{fit}:}{Object of class \code{"list"}, 
    	representing the MCMC draws. }
    \item{\code{summary}:}{Object of class \code{"list"}, 
    	representing the summary statistics. }
    \item{\code{setting}:}{Object of class \code{"list"}, 
    	representing the setting for graph-GPA 2.0 model fitting. }
    \item{\code{gwasPval}:}{Object of class \code{"matrix"}, 
    		representing the p-value matrix from GWAS data. }
    \item{\code{pgraph}:}{Object of class \code{"matrix"}, 
    		representing the prior phenotype graph. }
    \item{\code{annotMat}:}{Object of class \code{"matrix"}, 
    		representing the annotation matrix. }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GGPA2")}: provide brief summary of the object. }
    \item{plot}{\code{signature(x = "GGPA2", y = "missing", pCutoff = 0.5, betaCI = 0.95)}: plot a phenotype graph.
    Nodes i and j are connected if the posterior probability of E_{ij} > \code{pCutoff} and the posterior probability of beta_{ij} > \code{betaCI}. }
    \item{fdr}{\code{signature(object = "GGPA2", i=NULL, j=NULL)}: provide local FDR.
    By default (i.e., \code{i=NULL, j=NULL}),
    it returns a matrix of local FDR that a SNP is not associated with each phenotype (i.e., marginal FDR),
		where the order of columns is same as that in input GWAS data.
    If phenotype indices i and j are specified, a vector of corresponding local FDR is provided.
    }
    \item{estimates}{\code{signature(object = "GGPA2")}: 
        extract parameter estimates from graph-GPA 2.0 model fit. }
    \item{plotAnnot}{\code{signature(x = "GGPA2", namePheno = NULL, nameAnnot = NULL)}: plot boxplots of estimates of functional annotations. They are Markov random field coefficient estimates for gammas.}
     }

}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388
}
\author{ Hang J. Kim, Dongjun Chung, Qiaolan Deng }
\seealso{
\code{\link{GGPA2}}.
}
\examples{
showClass("GGPA2")

# Load the included simulation data
data(simulationAnnot)

# fit GGPA 2.0 model with 200 iterations and a burn-in of 200 iterations
# Note that we recommend more than 200 iterations in practice
fit <- GGPA2( gwasPval = simulationAnnot$pmat, annotMat = simulationAnnot$A_mat, nMain = 200, nBurnin = 200)

# Plot GGPA2 model fit
plot(fit, type="phenotype")

# Association mapping
head(fdr( fit ))
head(fdr( fit, i=1, j=2 ))
str(estimates( fit ))

# Plot boxplots of estimates of functional annotations
plot(fit, type="annotation")
}
\keyword{classes}
