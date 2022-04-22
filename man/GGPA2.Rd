\name{GGPA2}
\alias{GGPA2}
\title{
Fit graph-GPA 2.0 model
}
\description{
Fit graph-GPA 2.0 model.
}
\usage{
GGPA2( gwasPval, pgraph=NULL, annotMat = NULL, nBurnin=10000, nMain=40000, lbPval=1e-10, verbose=1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{gwasPval}{p-value matrix from GWAS data,
		where its row and column correspond to SNP and phenotype, respectively.}
    \item{pgraph}{A binary matrix representing the prior phenotype graph,
		where its rows and columns match the columns of \code{gwasPval}.}
    \item{annotMat}{A binary matrix representing the functional annotation data,
		where its row and column correspond to SNP and annotation, respectively.}
    \item{nBurnin}{Number of burn-in iterations. Default is 10000.}
    \item{nMain}{Number of main MCMC iterations. Default is 40000.}
    \item{lbPval}{Lower bound for GWAS p-value. 
  		Any GWAS p-values smaller than \code{lbPval} are set to \code{lbPval}.
  		Default is 1e-10.}  
    \item{verbose}{Amount of progress report during the fitting procedure. 
  		Possible values are 0 (minimal output), 1, 2, or 3 (maximal output).
  		Default is 1.} 
}
\details{
\code{GGPA2} fits the graph-GPA 2.0 model. It requires to provide GWAS p-value to \code{gwasPval}. If a phenotype graph is provided in \code{pgraph}, it is utilized to guide the phenotype graph estimation. Based on this GGPA2 fit, \code{assoc} implements association mapping and \code{plot} provides a phenotype graph.
}
\value{
Construct \code{GGPA2} class object.
}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388
}
\author{ Hang J. Kim, Dongjun Chung, Qiaolan Deng }
\seealso{
\code{\link{assoc}}, \code{\linkS4class{GGPA2}}.
}
\examples{


# Load the included simulation data without annotations
data(simulation)

# fit GGPA2 model with 200 iterations and a burn-in of 200 iterations
# Note that we recommend more than 200 iterations in practice
fit <- GGPA2( simulation$pmat, nMain = 200, nBurnin = 200)

# Association mapping with FDR of 0.1 and global control
head(assoc( fit, FDR=0.1, fdrControl="global" ))

# We may specift i = 1 and j = 2 if we are interested in that specific phenotype
head(assoc( fit, FDR=0.1, fdrControl="global", i=1, j=2 )) 

# plot the GGPA2 model fit
plot(fit)

# fit GGPA2 model with a prior graph
fit_prior <- GGPA2(simulation$pmat, pgraph = simulation$true_G, nMain = 200, nBurnin = 200)

# fit GGPA2 model with a prior graph and annotations
data(simulation_Annot)
fit_annot <- GGPA2(simulation_Annot$pmat, pgraph = simulation_Annot$true_G, annotMat = simulation_Annot$A_mat, nMain = 200, nBurnin = 200)


}
\keyword{models}
\keyword{methods}
