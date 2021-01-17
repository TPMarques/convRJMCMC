\name{plot.CZ_MANOVA}
\alias{plot.CZ_MANOVA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diagnostic graphics proposed by \insertCite{CZ;textual}{RJCA} for a single parameter}
\description{This function plots the 3 diagnostic graphics proposed by \insertCite{CZ;textual}{RJCA} for multiple parameters: 1 - MPSRF1 and MPSRF2 against batch; 2 - Largest eigenvalue of V and largest eigenvalue of Wc against batch; 3 - Largest eigenvalue of Wm and largest eigenvalue of WmWc against batch.}
\usage{\method{plot}{CZ_MANOVA}(x,...)}
%- maybe also 'usage' for other objects documented here.
\arguments{\item{x}{Object of class CZ_MANOVA created by CZ_MANOVA function.}\item{...}{Extra arguments to be parsed in the function.}}
\details{The plots are created in the currently active device. Users must press enter to see other plots.}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{\insertAllCited{}}
\author{Tiago Pereira Marques and Miriam Harumi Tsunemi}
\note{Axis and legend labels are in Brazilian Portuguese in this preliminary version. Future versions will have English labels as standard and allow user-defined labels.}

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

#CZ_MANOVA object is saved on test_2
teste_2<-CZ_MANOVA(theta,chains,models,mcmciterations)
plot(teste_2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{RJMCMC}
\keyword{Convergence}
\keyword{Graphics}
\keyword{Multiple parameters}
