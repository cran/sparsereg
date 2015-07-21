\name{volcanoplot}
\alias{volcanoplot}
\title{
Function for plotting posterior distribution of effects of interest.
}
\description{
The function produces a volcano plot for specified effects.  This can be useful for presenting or examining particular marginal effects of interest.
}
\usage{
volcanoplot(x, columns=NULL, newlabels=NULL, type="mode", stage=NULL)
}
\arguments{
  \item{x}{
Object of class \env{sparsereg}.
}
\item{columns}{ A vector of numbers (or strings) corresponding to columns (or column names) to produce plots for.}

\item{newlabels}{ New labels for columns rather than variable names in object.  If empty, variable names are used.}

\item{type}{ Options are \env{"mode"} and \env{"mean"}. Whether to plot the posterior mode or mean.}

\item{stage}{  When plotting output from a sparsereg object, the argument is ignored. When plotting output from a type 2 tobit, whether to return stage 1 or stage 2 results.}

}



\details{
Generates a volcano plot for coefficients from object from class \env{sparsereg}.  The desired coefficients can be requested using the \env{columns} argument and they can be assigned new names through \env{newlabels}.   }
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}
\seealso{\link{sparsereg}, \link{plot.sparsereg}, \link{summary.sparsereg}, \link{difference},  \link{print.sparsereg}}

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
 volcanoplot(s1,1:3)
}

}
