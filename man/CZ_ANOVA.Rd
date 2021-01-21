\name{CZ_ANOVA}
\alias{CZ_ANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RJMCMC convergence assessment of a single parameter
%%  ~~function to do ... ~~
}
\description{
This function creates an object of class CZ_ANOVA with the necessary information to plot the diagnostic graphics proposed by the \insertCite{CZ;textual}{convRJMCMC} criteria for a single parameter.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CZ_ANOVA(theta,chains,models,mcmciterations,nbatches,batchsize,confidence,division)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
  Observed values for the parameter in the different chains.
%%     ~~Describe \code{x} here~~
}
  \item{chains}{
  Chain index for the observed parameter value. The number of iterations per chain must be equal.
  }
  \item{models}{
  Model index for the observed parameter value.
  }
  \item{mcmciterations}{
  For a given chain is the iteration index value.
  }
  \item{nbatches}{
  Numbers of batches. The standard value is 20.
  }
  \item{batchsize}{
  It is the size of a batch given by iterations in a chain/nbatches consecutive integer. The function disregard nbatches value in case both are available. There is only a need to supply the nbatches or batchsize value.
  }
  \item{confidence}{
 It is the confidence level for the upper bound of both potential scale reduction factors. The standard value is 0.95.
  }
  \item{division}{
  Division criteria to make the batches. If division="Batch" the batches are partitions with size=batchsize of the iterations set. If division="Sequential" all the batches will have size=number of iterations/2, with the first batch starting at first iteration and the others at iteration 1+(number of iterations/2)/nbatches. Standard is division="Batch".
  }
}
\details{
The function output can be different from \insertCite{CZ;textual}{convRJMCMC} criteria in cases that at least one model is visited only for some chains as the degrees of freedom are automatically calculated by aov function. We recommend not to use the original criteria from \insertCite{CZ;textual}{convRJMCMC} in the situation above and to use this function instead.
}
\value{
An object of class CZ_ANOVA which saves a list for every batch with the summary of the ANOVA 1 model, the summary of the ANOVA 2 Model, the summary of the ANOVA 3 model, PSRF1, PSRF2, the upper bound of PSRF1, the upper bound of PSRF2, V, Wc, Wm, and WcWm. For more details see \insertCite{CZ;textual}{convRJMCMC} article.
}
\references{
  \insertAllCited{}
%% ~put references to the literature/web site here ~
}
\author{
Tiago Pereira Marques and Miriam Harumi Tsunemi
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#Example 1

set.seed(634963)
#Create parameters observations for the 3 chains
#theta 1 - chains 1 ...
theta1<-rnorm(100000,mean=20,sd=1.5)
theta2<-rnorm(100000,mean=20,sd=1.5)
theta3<-rnorm(100000,mean=20,sd=1.5)
#Merge generate observations in theta object
theta<-c(theta1,theta2,theta3)
#Create chain index
chains<-c(rep(1,100000),rep(2,100000),rep(3,100000))
#Create model index
  #In these configuration there will be 3 models in
  #proportions 1:2:1
models<-factor(sample(c(1,2,3),300000,TRUE,c(0.25,0.5,0.25)))
#Create iterations index
mcmciterations<-rep(seq(1:100000),3)
#Create CZ_ANOVA object and save to test 1
teste_1<-CZ_ANOVA(theta,chains,models,mcmciterations)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RJMCMC}% use one of  RShowDoc("KEYWORDS")
\keyword{Convergence}
\keyword{Single parameter}% __ONLY ONE__ keyword per line

