\name{makeSimul}
\alias{makeSimul}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate simulated data for graph-GPA 2.0 model
}
\description{
Generate simulated data for graph-GPA 2.0 model.
}
\usage{
makeSimul(true_G, Beta_ij, Alpha_i, mu1, sigma1, true_Gamma = NULL, A_mat = NULL, n_SNP=20000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{true_G}{
  An n by n matrix represents the true pairwise relationship between phenotypes (\eqn{E_{ij}}). 
  }
  \item{Beta_ij}{
  An n by n matrix represents the true \eqn{\beta_{ij}}. \eqn{\beta_{ij}} is the MRF coefficient for the pair of phenotypes i and j.
  }
  \item{Alpha_i}{
  A vector of length n represents the true \eqn{\alpha_{i}}. \eqn{\alpha_{i}} is the MRF coefficient for the phenotype i.
  }
  \item{mu1}{
  A vector of length n represents the true \eqn{\mu_{i}}. \eqn{\mu_{i}} is the mean of the log normal distribution.
  }
  \item{sigma1}{
  A vector of length n represents the true \eqn{\sigma_{i}}. \eqn{\sigma_{i}} is the variance of the log normal distribution.
  }
  \item{true_Gamma}{
  An n by M matrix represents the true \eqn{\gamma_{im}}. \eqn{\gamma_{im}} is the coefficient for importance of annotation m for phenotype i. If Null, simulated data without annotations would be generated. true_Gamma and A_mat must be set jointly.
  }
  \item{A_mat}{
  An M by T matrix consists of binary indicators representing the true annotations. 1 means the SNP is functional in the corresponding annotation. If Null, simulated data without annotations would be generated. true_Gamma and A_mat must be set jointly.
  }
  \item{n_SNP}{
  Number of SNPs to generate in the simulated dataset. Default is 20,000.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Construct a list of simulated data.
}
\references{
Deng Q, Nam JH, Yilmaz AS, Chang W, Pietrzak M, Li L, Kim HJ, and Chung D (2022), "graph-GPA 2.0: A graphical model for multi-disease analysis of GWAS results with integration of functional annotation data."

Kim H, Yu Z, Lawson A, Zhao H, and Chung D (2018), "Improving SNP prioritization and pleiotropic architecture estimation by incorporating prior knowledge using graph-GPA," *Bioinformatics*, 34(12): 2139-2141.

Chung D, Kim H, and Zhao H (2016), 
"graph-GPA: A graphical model for prioritizing GWAS results and investigating pleiotropic architecture," 13(2): e1005388
}
\author{
Hang J. Kim, Dongjun Chung, Qiaolan Deng 
}


\examples{
data("simulation")
dim(simulation$pmat)

# The simulated data is going to include 6 phenotypes.
# true beta and alpha
true_beta <- simulation$true_beta[1:6,1:6]
true_beta[1,2] <- true_beta[2,1] <- 3.2


# Generate the annotation matrix, A_annot, with 5 annotations.
n = 6
T = 20000
M = 5

#set.seed(1)
A_annot = matrix(rbinom(T*M,size=1,prob = 0.1),nrow = M)


# true mu and sigma
true_mu <- c(1.05,0.9,1.0,1.0,1.05,0.95)
true_sigma <- c(0.4,0.3,0.35,0.3,0.45,0.4)

# true Gamma
true_Gamma = rbind(
  c(1,0,0,0,0),
  c(1,0,0,0,0),
  c(1,0,0,0,0),
  c(0,2,0,0,0),
  c(0,2,0,0,0),
  c(0,2,0,0,0)
)

simdata <- makeSimul(true_G = adjmat, 
                  Beta_ij = true_beta, 
                  Alpha_i = diag(true_beta), 
                  mu1 = true_mu, 
                  sigma1 = true_sigma, 
                  true_Gamma = true_Gamma, 
                  A_mat = A_annot, 
                  n_SNP=T)
}
\keyword{methods}
