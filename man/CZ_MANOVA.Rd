\name{CZ_MANOVA}
\alias{CZ_MANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
RJMCMC convergence assessment of multiple parameters
}
\description{
This function creates an object of class CZ_MANOVA with the necessary information to plot the diagnostic graphics proposed by the \insertCite{CZ;textual}{RJCA} criteria for multiple parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CZ_MANOVA(theta,chains,models,mcmciterations,nbatches,batchsize,confidence,division)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
Matrix with parameters observations, different parameters must be on different collumns. For some iteration, chain, and model in a row, there must be a vector of the observed paramenters.
}
  \item{chains}{
  Chain index for the observed parameters. The number of iterations per chain must be equal.
%%     ~~Describe \code{chains} here~~
}
  \item{models}{
%%     ~~Describe \code{models} here~~
  Model index for the observed parameters.
}
  \item{mcmciterations}{
%%     ~~Describe \code{mcmciterations} here~~
  For a given chain is the iteration index value.
}
  \item{nbatches}{
%%     ~~Describe \code{nbatches} here~~
  Numbers of batches. The standard value is 20.
}
  \item{batchsize}{
%%     ~~Describe \code{batchsize} here~~
  It is the size of a batch given by iterations in a chain/nbatches consecutive integer. The function disregard nbatches value in case both are available. There is only a need to supply the nbatches or batchsize value.
}
  \item{confidence}{
%%     ~~Describe \code{confidence} here~~
  It is the confidence level for the upper bound of both potential scale reduction factors. The standard value is 0.95.
}
  \item{division}{
%%     ~~Describe \code{division} here~~
   Division criteria to make the batches. If division="Batch" the batches are partitions with size=batchsize of the iterations set. If division="Sequential" all the batches will have size=number of iterations/2, with the first batch starting at first iteration and the others at iteration 1+(number of iterations/2)/nbatches. Standard is division="Batch".
}
}
\details{
The function output can be different from \insertCite{CZ;textual}{RJCA} criteria in cases that at least one model is visited only for some chains as the degrees of freedom are automatically calculated by the manova function. We recommend not to use the original criteria from \insertCite{CZ;textual}{RJCA} in the situation above and to use this function instead.
}
\value{
An object of class CZ_MANOVA which saves a list for every batch with the V matrix, the Wc matrix, the Wm matrix, the WmWc matrix, MPSRF1, MPSRF2, the V matrix largest eigenvalue, the Wc matrix largest eigenvalue, the Wm matrix largest eigenvalue, and the WmWc matrix largest eigenvalue. For more details see \insertCite{CZ;textual}{RJCA} article.
}
\references{\insertAllCited{}}
\author{Tiago Pereira Marques and Miriam Harumi Tsunemi}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#Example 1

set.seed(500)

chains<-c(rep(1,100000),rep(2,100000),rep(3,100000))
models<-factor(sample(c(1,2,3),300000,TRUE,c(0.25,0.5,0.25)))
mcmciterations<-rep(seq(1:100000),3)

{
theta1a<-rnorm(100000,mean=0,sd=1)
theta2a<-rnorm(100000,mean=1,sd=2)
theta3a<-rnorm(100000,mean=1,sd=3)
thetaa<-cbind(theta1a,theta2a,theta3a)

theta1b<-rnorm(100000,mean=0,sd=1)
theta2b<-rnorm(100000,mean=1,sd=2)
theta3b<-rnorm(100000,mean=1,sd=3)
thetab<-cbind(theta1b,theta2b,theta3b)

theta1c<-rnorm(100000,mean=0,sd=1)
theta2c<-rnorm(100000,mean=1,sd=2)
theta3c<-rnorm(100000,mean=1,sd=3)
thetac<-cbind(theta1c,theta2c,theta3c)

theta<-rbind(thetaa,thetab,thetac)
}

teste_2<-CZ_MANOVA(theta,chains,models,mcmciterations)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RJMCMC}% use one of  RShowDoc("KEYWORDS")
\keyword{Convergence}
\keyword{Multiple parameters}

