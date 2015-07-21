# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

////// update beta //////

// [[Rcpp::export()]]
Rcpp::List updatebeta_cpp(Rcpp::NumericVector X0,
			  Rcpp::NumericVector y0,
			  Rcpp::NumericVector betacurr0,
			  Rcpp::NumericVector betamode0,
			  Rcpp::NumericVector lambdavec0,
			  Rcpp::NumericVector dtau0,
			  Rcpp::NumericVector sigmasq0,
			  Rcpp::NumericVector ps_sigmasq0,
			  Rcpp::NumericVector lambdashrink0,
			  Rcpp::NumericVector k0

			  ) {

//Declare inputs
arma::mat X = Rcpp::as< arma::mat >(X0);
arma::mat y = Rcpp::as< arma::mat >(y0);
arma::mat betacurr = Rcpp::as< arma::mat >(betacurr0);
arma::mat betamode = Rcpp::as< arma::mat >(betamode0);
arma::mat lambdavec = Rcpp::as< arma::mat >(lambdavec0);
arma::mat lambdashrink = Rcpp::as< arma::mat >(lambdashrink0);
arma::mat dtau = Rcpp::as< arma::mat >(dtau0);
arma::mat betaci = betamode;
arma::mat betals = betamode;

arma::mat sigmasq=Rcpp::as< arma::mat >(sigmasq0);
arma::mat ps_sigmasq=Rcpp::as< arma::mat >(ps_sigmasq0);
arma::mat k=Rcpp::as< arma::mat >(k0);


int p = X.n_cols;
int n = X.n_rows;

double delbeta=0;
arma::mat fits = X*betacurr;
arma::mat var_beta = betacurr;

//Start loop
for(int i_coef = 0; i_coef<p; i_coef++){

//fits = X*betacurr;
//Generate components
arma::mat x_curr = X.submat(0,i_coef,n-1,i_coef);
arma::mat v_curr = y-fits + x_curr*betacurr(i_coef,0);
//arma::mat A_curr=arma::strans(x_curr)*x_curr + 1/dtau(i_coef,0);
double A_curr= (n-1) + 1/dtau(i_coef,0);
arma::mat olsnum = arma::strans(x_curr)*v_curr;
//arma::mat olsdenom = arma::strans(x_curr)*x_curr;
double olsdenom = n-1;
double ols_est = olsnum(0,0)/olsdenom;
betals(i_coef,0) = ols_est;

//Update posterior mean and CI
double meantemp = olsnum(0,0)/A_curr;
//Rcpp::NumericVector ranscalar = rnorm(1,0,1);

//Update CI
//betaci(i_coef,0) =  (pow((1+p)*sigmasq(0,0)/olsdenom,.5)); 
 arma::mat errs = y-fits;
 arma::mat temp = arma::strans(errs)*errs;
 betaci(i_coef,0)=temp(0,0)/(n-1);

double oldbeta = betacurr(i_coef,0);
//double newbeta =meantemp+ranscalar(0)*pow(sigmasq(0,0)/A_curr,.5); 

var_beta(i_coef,0)= sigmasq(0,0)/A_curr;

//Update with antithetical gibbs sampler
double newbeta = 2*meantemp-oldbeta;

//Update with gibbs sampler
//double newbeta =meantemp+ranscalar(0)*pow(var_beta(i_coef,0),.5); 


betacurr(i_coef,0) =newbeta;
double delbeta = newbeta-oldbeta;

	for(int i_fit =0; i_fit<n; i_fit ++) 
			{
			fits(i_fit,0) = fits(i_fit,0) + delbeta*x_curr(i_fit);
			}


double ols_est_0=ols_est;

//Update posterior mode
if(ols_est>0) {
	//Bayesian LASSO
	ols_est = ols_est - k(0,0)*lambdashrink(i_coef,0)*pow(ps_sigmasq(0,0),.5)/pow(olsdenom,1);		
	if(ols_est <=0) ols_est=0;
	if(ols_est > 0) ols_est = betacurr(i_coef,0);


	//Generalized double exponential
	//ols_est = ols_est - 2*sigmasq(0,0)/(pow(sigmasq(0,0),.5)+abstemp )/pow(olsdenom,1);
	//ols_est = ols_est - 2*pow(ps_sigmasq(0,0),-.5)/((1+pow(ps_sigmasq(0,0),-.5)*abstemp))*pow(olsdenom,-1);
	//if(ols_est<0) ols_est=0;
		
	//double a = 0;
	//double b = 0;
	//double c = 0;
	
	//a = olsdenom/pow(ps_sigmasq(0,0),.5);
	//b = -olsdenom+olsnum(0,0)*pow(ps_sigmasq(0,0),-.5);
	//c = 2*pow(ps_sigmasq(0,0),-.5)-olsnum(0,0);

	//ols_est =  (-b+pow(b*b-4*a*c,.5))/(2*a);
	//ols_est = -ols_est;

	//if(b*b-4*a*c<0) ols_est=0;
	if(ols_est<0) ols_est = 0;
	if(ols_est!=ols_est) ols_est=0;
	if(betacurr(i_coef,0)<ols_est) ols_est=0;

	}
if(ols_est<0) {
	//Bayesian LASSO
	ols_est = ols_est + k(0,0)*lambdashrink(i_coef,0)*pow(ps_sigmasq(0,0),.5)/pow(olsdenom,1);	
	if(ols_est >=0) ols_est=0;
	if(ols_est < 0) ols_est = betacurr(i_coef,0);

	//Generalized double exponential
	//ols_est = ols_est + 2*pow(sigmasq(0,0),.5)/((1+pow(sigmasq(0,0),-.5)*abstemp)*pow(olsdenom,1));
	//if(ols_est_0 + 2*sigmasq(0,0)/(pow(sigmasq(0,0),.5)+abstemp )/pow(olsdenom,1) >=0) ols_est=0;
	//ols_est = ols_est + 2*pow(ps_sigmasq(0,0),-.5)/((1+pow(ps_sigmasq(0,0),-.5)*abstemp))*pow(olsdenom,-1);
	//if(ols_est>0) ols_est=0;

	//double a = 0;
	//double b = 0;
	//double c = 0;
	
	//a = -olsdenom/pow(ps_sigmasq(0,0),.5);
	//b = -olsdenom-olsnum(0,0)*pow(ps_sigmasq(0,0),-.5);
	//c = -2*pow(ps_sigmasq(0,0),-.5)-olsnum(0,0);

	//ols_est =  (-b-pow(b*b-4*a*c,.5))/(2*a);

	//ols_est = -ols_est;

	if(ols_est>0) ols_est =0;
	if(ols_est!=ols_est) ols_est=0;
	if(betacurr(i_coef,0)>ols_est) ols_est=0;

	}

int signcheck=1;
int signest=1;
if(ols_est<0) signest = -1;
if(betacurr(i_coef,0)<0) signcheck=-1;
if(signcheck* signest< 0) ols_est=0;
betamode(i_coef) = ols_est;
}		




return Rcpp::List::create(Rcpp::Named("beta.mean") = betacurr, Rcpp::Named("beta.mode") = betamode, Rcpp::Named("var.beta") = var_beta,Rcpp::Named("beta.ci")= betaci,
Rcpp::Named("beta.ols") = betals
);


}
