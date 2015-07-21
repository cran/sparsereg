\name{print.sparsereg}
\alias{print.sparsereg}
\title{
A summary of the estimated posterior mode of each parameter.
}
\description{
The funciton prints a summary of the estimated posterior mode of each parameter.
}
\usage{
\method{print}{sparsereg}(x,... )
}
\arguments{
  \item{x}{
Object of class \env{sparsereg}.
}
\item{...}{
	Additional arguments to pass to \env{print}.  None supported in this version.
	}
}


\details{Uses the \env{summary} function from the package \env{coda} to return a summary of the posterior mode of a \env{sparsereg} object.  The argument \env{stage} is required when printing output from a type 2 tobit.	 When printing output from a sparsereg object, the argument is ignored. 
}

	
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}
\seealso{\link{sparsereg}, \link{plot.sparsereg}, \link{summary.sparsereg}, \link{volcanoplot}, \link{difference}}
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
 print(s1)
}

}
