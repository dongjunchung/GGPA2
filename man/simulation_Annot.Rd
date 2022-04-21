\name{simulation_Annot}
\alias{simulation_Annot}
\docType{data}
\title{Simulation data for graph-GPA 2.0}
\description{
This is an simulation dataset with annotations.
}
\usage{data("simulation_Annot")}
\format{
  \code{simulation_Annot} list object containing simulation data (element Y) and its simulation setting (the remaining elements).
}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388
}
\author{ Hang J. Kim, Dongjun Chung, Qiaolan Deng }
\examples{

# The simulation data set is included with the GGPA2 package
data(simulation_Annot)
head(t(simulation_Annot$pmat))
head(t(simulation_Annot$A_mat))
}
\keyword{datasets}
