\name{assoc}
\alias{assoc}
\alias{assoc,GGPA2-method}
\title{
Association mapping
}
\description{
Association mapping.
}
\usage{
assoc( object, ... )
\S4method{assoc}{GGPA2}( object, FDR=0.05, fdrControl="global", i=NULL, j=NULL )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{object}{A GGPA2 model fit as obtained by GGPA2(). }
    \item{FDR}{The desired FDR level. }
	  \item{fdrControl}{Method to control FDR. Possible values are "global" (global FDR control) and
        "local" (local FDR control). Default is "global". } 
    \item{i}{Index for the first phenotype used in association mapping.
      See the details about how users can specify the pattern.}
    \item{j}{Index for the second phenotype used in association mapping.
      See the details about how users can specify the pattern.}
    \item{...}{ Other parameters to be passed through to generic \code{assoc}.}
}
\details{
\code{assoc} uses the direct posterior probability approach of Newton et al. (2004)
to control global FDR in association mapping.

By default (i.e., \code{i=NULL, j=NULL}), \code{assoc} implements association mapping for each phenotype. If users are interested in identifying SNPs associated with a pair of phenotypes, users can specify indices of phenotypes of interest using the arguments \code{i} and \code{j}. Note that both \code{i} and \code{j} should be either NULL or numeric.
}
\value{
If \code{i=NULL, j=NULL}, returns a binary matrix indicating association of SNPs for each phenotype,
where its rows and columns match those of input p-value matrix for function \code{GGPA}.
Otherwise, returns a binary vector indicating association of SNPs for i-th and j-th phenotype pair.
}
\author{
Hang J. Kim, Dongjun Chung, Qiaolan Deng
}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388

Newton MA, Noueiry A, Sarkar D, and Ahlquist P (2004),
"Detecting differential gene expression with a semiparametric hierarchical mixture method,"
\emph{Biostatistics}, Vol. 5, pp. 155-176.
}
\seealso{
\code{\link{GGPA2}}, 
\code{\linkS4class{GGPA2}}.
}
\examples{

# Load the included simulation data
data(simulation)

# fit GGPA2 model with 200 iterations and a burn-in of 200 iterations
# Note that we recommend more than 200 iterations in practice
fit <- GGPA2( simulation$pmat, nMain = 200, nBurnin = 200)

# Association mapping with FDR of 0.1 and global control
head(assoc( fit, FDR=0.1, fdrControl="global" ))

# We may specift i = 1 and j = 2 if we are interested in that specific phenotype
head(assoc( fit, FDR=0.1, fdrControl="global", i=1, j=2 )) 

}
\keyword{models}
\keyword{methods}
