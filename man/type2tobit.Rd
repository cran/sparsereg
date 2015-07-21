\name{type2tobit}
\alias{type2tobit}
\title{Type 2 Tobit model using sparse regression.}
\description{Function for fitting a type 2 Tobit model.  Both stages of the model are fit using sparse modeling.}
\usage{
type2tobit(y1, y2=NULL, X1, X2=NULL, stage2, 
group=NULL, gibbs=200, burnin=200, thin=10,
one.constraint=TRUE, id=NULL, id2=NULL, id3=NULL, 
id_2=NULL, id2_2=NULL, id3_2=NULL)
}
\arguments{
  \item{y1, y2}{Dependent variables for stage 1 (binary model) and stage 2 (truncated model).}
\item{X1, X2}{Covariates for stage 1 (binary model) and stage 2 (truncated model).}
\item{stage2}{An indicator variable of length \env{y2} taking a value of 1 for observations for which the inverse mills ratio in the stage 2 model is zero.}
\item{group}{Optional vector of integers characterizing which covariates are under the same LASSO constraint.  For example, the vector c(1,1,2,2,2) will place the first two covariates under one LASSO constraint and the last three under another.}
\item{gibbs}{Number of posterior samples to save.  Between each saved sample,  \env{thin} samples are drawn.}
\item{burnin}{Number of burnin samples.   Between each burnin sample,  \env{thin} samples are drawn.
These iterations will not be included in the resulting analysis.}
\item{thin}{Extent of thinning of the MCMC chain.  Between each posterior sample, whether burnin or saved, \env{thin} draws are made.}
\item{one.constraint}{Whether to fit a single LASSO constraint.}
\item{id, id2, id3, id_2, id2_2, id3_2}{Vectors the same lenght of the sample denoting clustering in the data.  In a conjoint experiment with repeated observations, these correspond with respondent IDs.  Up to three different sets of random effects are allowed.}
}
\details{
Runs a type 2 tobit model with stage 1 running a probit and stage 2 entering the residual correction as a covariate.  Please note that if there are any known cases that have non-zero value with probability 1, these observations \emph{must} be stacked on top of the stage 1 covariates when they are given to stage 2.  If there are no such obserations, then stage 1 and stage 2 should consist of the same observations in the same order.
}
\value{
  \item{beta.mode}{Matrix of sparse (mode) estimates for stage 1 with columns equal to number of effects and rows for posterior samples.}
  \item{beta.mean}{Matrix of mean estimates for stage 1 with columns equal to number of effects and rows for posterior samples.  These estimates are not sparse, but they do predict better than the mode.}
  \item{beta.mode.2}{Matrix of sparse (mode) estimates for stage 2 with columns equal to number of effects and rows for posterior samples.}
  \item{beta.mean.2}{Matrix of mean estimates for stage 2 with columns equal to number of effects and rows for posterior samples.  These estimates are not sparse, but they do predict better than the mode.}
  \item{sigma.sq_2}{Vector of posterior estimate of error variance for stage 2.}
  \item{X}{Matrix of covariates fit.  Includes interaction terms, depending on \env{scale.type}.}
  \item{modeltype}{Type of sparsereg model fit. In this case, \env{twostage}.  Used by summary functions.}

 }
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.

Kim, In Song, Londregan, John and Ratkovic, Marc.  2015.  "A Sparse Decomposition of Democracy and Trade."  Working paper
}
\seealso{\link{sparsereg}}
\examples{
	\dontrun{
		
set.seed(1)
 n<-1000
 k<-5
 Sigma<-diag(k)
 Sigma[Sigma==0]<-0.5
 X<-mvrnorm(n,m=rep(0,k),S=Sigma)

 y1.0<- -1+X[,2]*1 + X[,3]*(-1)+rnorm(n)
 probs1<-dnorm(y1.0)
 y1<-obs.zero<-rbinom(n,1,probs1)

 y2<-exp(1+X[,1]*1 + X[,2]*(-2)+X[,4]*(-1)+rnorm(n))
 y2[obs.zero==0]<-0

 t1<-type2tobit(y1=y1,X1=X,y2=log(1+y2),X2=X,stage2=rep(0,n),
   one.constraint=TRUE)

 plot(t1, stage=1)
 plot(t1, stage=2)
		
	}
	
}

