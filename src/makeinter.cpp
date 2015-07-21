# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

////// make interactions: interact X0 and X1, place in X2 //////

// [[Rcpp::export()]]
Rcpp::List makeinter_cpp(Rcpp::NumericVector X00,
			  Rcpp::NumericVector X10,
  			  Rcpp::NumericVector X20
			  ) {

//Declare inputs
arma::mat X0 = Rcpp::as< arma::mat >(X00);
arma::mat X1 = Rcpp::as< arma::mat >(X10);
arma::mat X2 = Rcpp::as< arma::mat >(X20);


int p1 = X0.n_cols;
int p2 = X1.n_cols;
int n = X1.n_rows;

arma::mat Xtemp=arma::zeros<arma::mat>(n,3);
arma::mat XtX = arma::strans(Xtemp)*Xtemp;



for(int idum = 0; idum < n; idum++) Xtemp(idum,0) =1;

int pcurr = 0;
 
for(int i1 = 0; i1<p1; i1++){
for(int i2 = 0; i2<p2; i2++){

//Generate components
arma::mat x0_curr = X0.submat(0,i1,n-1,i1);
arma::mat x1_curr = X1.submat(0,i2,n-1,i2);
arma::mat inter = x1_curr;

//Generate interaction term
for(int idum = 0; idum < n; idum++) {
	inter(idum,0) =x0_curr(idum,0)*x1_curr(idum,0);
	Xtemp(idum,1)=x0_curr(idum,0);
	Xtemp(idum,2)=x1_curr(idum,0);
	}

//Generate hat matrix

XtX = arma::strans(Xtemp)*Xtemp;
XtX=.5*XtX+.5*arma::strans(XtX);
arma::mat XtXinvXty=arma::pinv(XtX)*(arma::strans(Xtemp)*inter);
inter = inter-Xtemp*XtXinvXty;

//Place in X2
for(int idum = 0; idum < n; idum++) {
	X2(idum,pcurr) =inter(idum,0);
	}

pcurr = pcurr+1;

}	
	
}


return Rcpp::List::create(Rcpp::Named("X") = X2);


}
