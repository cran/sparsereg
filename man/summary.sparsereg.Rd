\name{summary.sparsereg}
\alias{summary.sparsereg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summaries for a sparse regression.
}
\description{
	The function prints and returns a summary table for a \env{sparsereg} object.
	
}
\usage{
\method{summary}{sparsereg}(object,... )
}
\arguments{
  \item{object}{
Object of type \env{sparsereg}.
}
\item{...}{Additional items to pass to \env{summary}.  Options below.}

}



\details{
Generates a table for an object of class \env{sparsereg}.  Additional arguments to pass \env{summary} below.

\env{interval}{
	Length of posterior interval to return.  Must be between 0 and 1, default is .9.  The symmetric interval is returned.
}

\env{ci}{
	Type of interval to return.  Options are \env{"quantile"} (default) for quantiles and \env{"HPD"} for the highest posterior density interval.
}

\env{order}{
How to order returned coefficients.  Options are \env{"magnitude"}, sorted by magnitude and omitting zero effects, \env{"sort"}, sorted by size from highest to lowest and omitting zero effects,  and \env{"none"} which returns all effects
}

\env{normal}{
Whether to return the normal approximate confidence interval (default of TRUE) or posterior interval (FALSE).
}

\env{select}{
	Either \env{"mode"} or a number between 0 and 1.  Whether to select variables for printing off the median of the mode (default) or off the probability of being non-zero.  
}

\env{printit}{
Whether to print a summary table.}

\env{stage}{
	 When summarizing output from a sparsereg object, the argument is ignored. When summarizing output from a type 2 tobit, whether to return stage 1 or stage 2 results.
}

}
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}
\seealso{\link{sparsereg}, \link{plot.sparsereg}, \link{volcanoplot}, \link{difference}, \link{print.sparsereg}}
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
 summary(s1)
}

}
