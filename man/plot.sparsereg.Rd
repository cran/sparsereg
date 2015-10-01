\name{plot.sparsereg}
\alias{plot.sparsereg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting output from a sparse regression.
}
\description{
Function for plotting coefficients from sparsereg analysis.
}
\usage{
\method{plot}{sparsereg}(x,...)
}
\arguments{
  \item{x}{
Object from output of class \env{sparsereg}.
}
\item{...}{Additional items to pass to \env{plot}.  Options below.}
}
\details{
The function returns up to three plots in one figure.  Each plot corresponds with main effects, interaction effects, and two-way interactions.  Additional options to pass below.

\env{main1, main2, main3}{
Main titles for plots of main effects, interactive effects, and two-way interactions.}

\env{xlabel}{
Label for x-axis.
}

\env{plot.one}{
Takes on the value of FALSE or 1, 2, or 3, denoting whether to return a single plot for main effects (1), interactive effects (2), or two-way interactions (3).
}


}
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}

\seealso{\link{sparsereg}, \link{plot.sparsereg}, \link{summary.sparsereg}, \link{violinplot}, \link{difference}}
\examples{
\dontrun{
 set.seed(1)
 n<-500
 k<-100
 Sigma<-diag(k)
 Sigma[Sigma==0]<-.5
 X<-mvrnorm(n,mu=rep(0,k),Sigma=Sigma)
 y.true<-3+X[,2]*2+X[,3]*(-3)
 y<-y.true+rnorm(n)



##Fit a linear model with five covariates.
 s1<-sparsereg(y,X[,1:5])
 plot(s1)
}
}
