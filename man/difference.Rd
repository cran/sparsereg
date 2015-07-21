\name{difference}
\alias{difference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	
	Plotting difference in posterior estimates from a sparse regression.
}
\description{
Function for plotting differences in posterior density estimates for separate parameters from sparse regression analysis.
}
\usage{
difference(x,type="mode",var1=NULL,var2=NULL,plot.it=TRUE, 
main="Difference",xlabel="Effect", ylabel="Density")
}
\arguments{
  \item{x}{
Object of class \env{sparsereg}.
}
\item{type}{
	Whether to difference the posterior mode or posterior mean.  Options are \env{"mode"} and \env{"mean"}.
	}

\item{var1, var2}{
Variables names for the effects to difference.}
\item{plot.it}{Whether to plot the density of the difference.}
\item{main, xlabel, ylabel}{Main title, x-axis label, and y-axis label.}
	
}

\details{Generates a density of the estimated posterior of the difference between the effects of two variables. }
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}
\seealso{\link{sparsereg}, \link{plot.sparsereg}, \link{summary.sparsereg}, \link{volcanoplot},  \link{print.sparsereg}}
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
 difference(s1,var1=1,var2=2)
}

}
