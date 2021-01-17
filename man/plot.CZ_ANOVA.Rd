\name{plot.CZ_ANOVA}
\alias{plot.CZ_ANOVA}
\title{Diagnostic graphics proposed by \insertCite{CZ;textual}{RJCA} for a single parameter}
\description{This function plots the 3 diagnostic graphics proposed by \insertCite{CZ;textual}{RJCA} for a single parameter: 1 - PSRF1 and PSRF2 plus their upper bounds against batch; 2 - V and Wc against batch; 3 - Wm and WmWc against batch.}
\usage{\method{plot}{CZ_ANOVA}(x,...)}
\arguments{\item{x}{Object of class CZ_ANOVA created by CZ_ANOVA function.}\item{...}{Extra arguments to be parsed in the function.}}
\details{
The plots are created in the currently active device. Users must press enter to see other plots.}
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

set.seed(634963)
theta1<-rnorm(100000,mean=20,sd=1.5)
theta2<-rnorm(100000,mean=20,sd=1.5)
theta3<-rnorm(100000,mean=20,sd=1.5)
theta<-c(theta1,theta2,theta3)
chains<-c(rep(1,100000),rep(2,100000),rep(3,100000))
models<-factor(sample(c(1,2,3),300000,TRUE,c(0.25,0.5,0.25)))
mcmciterations<-rep(seq(1:100000),3)

#CZ_ANOVA object is saved on test_1
teste_1<-CZ_ANOVA(theta,chains,models,mcmciterations)
plot(teste_1)
}
\keyword{RJMCMC}
\keyword{Convergence}
\keyword{Graphics}
\keyword{Single parameter}

