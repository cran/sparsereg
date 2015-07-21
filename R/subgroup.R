
########################################
########################################
####The functions
########################################
########################################

#Takes type linear, tobit, probit

#scale.type: none, expand.treat, expand.X

sparsereg<-function(y, X, treat=NULL, gibbs=200, burnin=200, thin=10,type="linear", group=NULL,  trunc=NULL, lower.trunc=NULL, one.constraint=FALSE, scale.type="none", baseline.vec=NULL, id=NULL, id2=NULL, id3=NULL,save.temp=FALSE){

	##Warnings
	if(!scale.type%in%c("none","TX","TTX","TT")) stop("scale.type should be one of none, TX, TTX, or TT")
	if(!type%in%c("linear","probit","tobit")) stop("type should be one of linear, probit, or tobit")
	if(sum(!is.finite(y))>0) stop("Please remove missing or infinite values from y")
	if(sum(!is.finite(X))>0) stop("Please remove missing or infinite values from X")
	#if(nrow(as.matrix(X))!=length(y)) stop("y and X need to have the same length")
	if(length(id)>0& length(id)!=length(y)) stop("id and y need to have the same length")
	if(length(id2)>0& length(id2)!=length(y)) stop("id2 and y need to have the same length")
	if(length(id3)>0& length(id3)!=length(y)) stop("id3 and y need to have the same length")

	###################################
	###################################
	#### Declare variables
	###################################
	###################################
	#set.group<-length(group)!=0
	
	#if(length(X)>0) if(length(group)==0) group<-ifelse(length(X)>0,rep(1,ncol(X)),999)
	#if(length(group)!=0) group<-as.numeric(as.factor(group))
	
	id<-as.numeric(as.factor(id))
	id2<-as.numeric(as.factor(id2))
	id3<-as.numeric(as.factor(id3))

	if(length(treat)!=0) treat<-data.frame(treat)
	n<-length(y)

cat("Step 1 of 4: Formatting Data","\n")

	if(length(X)>0){
		if(length(colnames(X))==0) colnames(X)<-paste("X",1:ncol(X),sep="")
	X<-X[,apply(X,2,sd)>0]
	}
	
	if(type=="probit") {
	y0<-(y>mean(y))*1
	y[y0==1]<-dnorm(0)/pnorm(0)
	y[y0==0]<--dnorm(0)/pnorm(0)
	}
	
	y<-as.vector(y)
	if(length(X)>0){
		X<-as.matrix(X)
		X<-apply(X,2,as.numeric)
		}
	if(length(treat)==0) X.big<-X.lin<-X
	drop.string<-"IPREFERONEDIRECTIONTOLEDZEPPELIN"
	if(length(treat)!=0){
		treat<-as.matrix(treat)
		if(length(baseline.vec)>0){
		for(i.b in 1:length(baseline.vec)) 
		treat[treat[,i.b]==baseline.vec[i.b],i.b]<-drop.string
		}
			
	dummy.col<-function(t1){
		t1<-as.matrix(t1)
		if(length(colnames(t1))==0) colnames(t1)<-"treat"
		out<-sapply(sort(unique(t1)),FUN=function(x) x==t1)*1
		colnames(out)<-paste(colnames(t1)[1],sort(unique(t1)), sep="_")
		invisible(out)
	}
	if(length(treat)>0) {
	treat0<-as.matrix(treat)
	treat.mat.0<-NULL
	for(i.t in 1:ncol(treat)) treat.mat.0<-cbind(treat.mat.0,dummy.col(treat0[,i.t]))
	}
	if(length(X)>0) X.lin<-cbind(X,treat.mat.0) else X.lin<-treat.mat.0

	X.lin<-apply(X.lin,2,as.numeric)
	colnames(X.lin)<-gsub(":","_",colnames(X.lin))
	}
	drop.cols<-grep(drop.string,colnames(X.lin))
	if(length(drop.cols)>0) X.lin<-X.lin[,-drop.cols]
	
	
if(scale.type=="none"){
	if(length(X)>0) c.cluster<-c(rep(1,ncol(X)),rep(2,ncol(X.lin)-ncol(X))) else
		c.cluster<-rep(1,ncol(X.lin))
	X.lin<-apply(X.lin,2,FUN=function(x) x-mean(x) )
	X.big<-X.lin
	if(one.constraint) c.cluster<-rep(1,ncol(X.big))
	}


	if(scale.type=="TTX"){
			treat<-as.matrix(treat)
			all.vars<-make.threewayinter(X,treat,treat)
			X.big<- all.vars$big.X
			c.cluster<-all.vars$c.clust
	}


	if(scale.type=="TX"){
			if(length(X)==0) X<-matrix(1,nrow=length(y),ncol=1)
			treat<-as.matrix(treat)
			all.vars<-make.inter(X,treat)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust

	}

	if(scale.type=="TT"&(length(X)==0)){
			treat<-as.matrix(treat)
			X.lin<-apply(X.lin,2,as.numeric)
			all.vars<-make.inter(X.lin,X.lin)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust

	}


	if(scale.type=="TT"&(length(X)>0)){
			all.vars<-make.inter.treat(X,treat)
			X.big<- all.vars$X
			c.cluster<-all.vars$c.clust
		}


	drop.cols<-grep(drop.string,colnames(X.big))
	if(length(drop.cols)>0) X.big<-X.big[,-drop.cols]
	
	###################################
	###################################
	#### Eliminate correlated vars from X
	###################################
	###################################
	X.big[abs(X.big)<1e-10]<-0

	X.big<-apply(X.big,2,FUN=function(x) x-mean(x))
	keeps<-apply(X.big,2,sd)>0&(apply(X.big,2,FUN=function(x) sum(is.na(x)))==0)

	#if(set.group)	{c.cluster<-group
	#if(length(group)!=ncol(X.big)) stop("group should have one entry for every element of the full X matrix")
	#}

	X.big<-X.big[,keeps]
	c.cluster<-c.cluster[keeps]

	
	cor1<-cor(X.big)
	diag(cor1)<-0
	cors.run<-NULL
for(i in 1:(ncol(X.big)-1)) for(j in (i+1):ncol(X.big)){
	if(cor1[i,j]^2>0.999) cors.run<-rbind(cors.run,c(i,j))
}

	if(length(cors.run)>0){
	drops<-as.vector(cors.run[,2])
	X.big<-X.big[,-drops]
	c.cluster<-c.cluster[-drops]
	}
	

	scale.big<-sapply(colnames(X.lin),grep,colnames(X.big))
	matrix.convert<-matrix(0,nrow=length(colnames(X.lin)),ncol=ncol(X.big))
	rownames(matrix.convert)<-colnames(X.lin)
	colnames(matrix.convert)<-colnames(X.big)
	for(i.scale in 1:ncol(X.lin)) matrix.convert[i.scale,unlist(scale.big[i.scale])]<-1

	scale.back<-apply(X.big,2,FUN=function(x) sd(x))

	for(i.scale in 1:ncol(X.big)) X.big[,i.scale]<-X.big[,i.scale]/scale.back[i.scale]
	
	scale.y<-sd(y)
	y<-(y-mean(y))/scale.y
	if(type=="probit") {
	scale.y<-1
	y0<-(y>mean(y))*1
	y[y0==1]<-dnorm(0)/pnorm(0)
	y[y0==0]<--dnorm(0)/pnorm(0)
	}

	if(type=="probit"|type=="tobit"){
		X.big<-cbind(1,X.big)
		colnames(X.big)[1]<-"intercept"
		c.cluster<-c(max(c.cluster)+1,c.cluster)
		scale.back<-c(1,scale.back)
	}
	n.cluster<-length(unique(c.cluster))


	###################################
	###################################
	#### Format containers; initialize
	###################################
	###################################

	n<-length(y)
	p<-ncol(X.big)
	c.cluster<-c.cluster[1:p]

	beta.mean<-matrix(NA,nrow=p,ncol=gibbs)
	beta.mode<-matrix(NA,nrow=p,ncol=gibbs)
	lambda.run<-matrix(NA,nrow=p,ncol=gibbs)
	beta.ci<-matrix(NA,nrow=p,ncol=gibbs)

	sigma.sq.run<-rep(NA,gibbs)
	rownames(beta.mean)<-rownames(beta.mode)<-colnames(X.big)
	beta.curr.mode<-beta.curr<-rnorm(p,sd=.1)
	

	ps.sigmasq<-sigma.sq<-mean((y-X.big%*%beta.curr)^2)
	if(type=="probit") sigma.sq<-1

	lambda.shrink<-lambda.all<-rep(0,p)
	for(i in unique(c.cluster)) {
		lambda.all[c.cluster==i]<-(2*sum(c.cluster==i)/sum(abs(beta.curr)[c.cluster==i]))^.5
		}
	lambda.shrink<-lambda.all
	if(one.constraint) lambda.all[1:p]<-(2*p/sum(abs(beta.curr)))^.5
	k.vec<-rep(NA,length(c.cluster))
	for(i in 1:n.cluster) k.vec[c.cluster==i]<-sum(c.cluster==i)
	#if(ncol(X.big)<150) block<-TRUE
	#block<-FALSE
	#if(block) XprimeX<-t(X.big)%*%X.big


	##Initialize Random Effects
	ran.effect<-rep(0,n)
	sigma.sq.b<-1
	k.raneff<- length(unique(id))+length(unique(id2))+length(unique(id3))

	if(k.raneff>0){
	ran.effect1<-ran.effect2<-ran.effect3<-rep(0,n)
	sigma.sq.b<-sigma.sq.b2<-sigma.sq.b3<-1

	re1<-updateREs(y, fits.main=as.vector(X.big%*%beta.curr), ran.effect1, ran.effect2, ran.effect3, id, id2, id3,sigma.sq.b, sigma.sq.b2,sigma.sq.b3,sigma.sq=1,fix.eff=FALSE)
	sigma.sq.b<-re1$sigmas[1]
	sigma.sq.b2<-re1$sigmas[2]
	sigma.sq.b3<-re1$sigmas[3]
	ran.effect1<-re1$REs[,1]
	ran.effect2<-re1$REs[,2]
	ran.effect3<-re1$REs[,3]
	ran.effect<-rowSums(re1$REs)
	} 

cat("Step 2 of 4: Beginning Burnin Period","\n")

	##Begin gibbs loop
	for(i.gibbs in 1:(burnin+gibbs)){	
	for(i.thin in 1:thin){#Start thin loop
	##Update tau^2; need these for betas
	muprime.all<-abs(lambda.all*sqrt(sigma.sq))/abs(beta.curr)
	lambda.prime<-lambda.all^2
	invTau2<-sapply(1:length(muprime.all), FUN=function(i) rinv.gaussian(1, muprime.all[i], (lambda.prime[i] ) ) )
	Dtau<-abs(1/invTau2)
	
	up1<-updatebeta_cpp(X0=X.big, y0=as.matrix(y-ran.effect), betacurr0=as.matrix(beta.curr), betamode0=as.matrix(beta.curr.mode), lambdavec0=as.matrix(lambda.all), dtau0=as.matrix(Dtau), sigmasq0=as.matrix(sigma.sq),ps_sigmasq0=as.matrix(ps.sigmasq),lambdashrink0=as.matrix(lambda.shrink),
	k0=as.matrix(1)
	)
	beta.curr<-up1$beta.mean
	beta.curr.mode<-up1$beta.mode
	beta.curr.ci<-up1$beta.ci
	if(i.thin==thin){

#	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)

	for(i.var in 1:p){
        beta.ls<-up1$beta.ols[i.var]
        lambda.use<-lambda.all[i.var]
	in.phi<-abs(n^.5*abs(beta.ls/sigma.sq^.5)-lambda.shrink[i.var]*ps.sigmasq^.5/(sigma.sq^.5*(n-1)) ) 
	p.z<-pnorm(in.phi)
        beta.2<-beta.curr.mode[i.var]
	var.beta2<-sigma.sq/(n-1)
        var.betals<-beta.curr.ci[i.var]
	var.lambda<-1/(4*lambda.use^2)*(p*n^.5+1)/(sum(Dtau)/2+1.78)^2
        beta.sigma<-sigma.sq.scale#(beta.curr.ci[i.var]*(n-1)+sum(beta.curr^2/Dtau))/2
        alpha.sigma<-ps.sigma.sq.shape#(n^.5-1)/2
        #Variance of pseudo sigma
        var.sigma<-1/(4*ps.sigmasq)*beta.sigma/((alpha.sigma-1)^2*(alpha.sigma-2))
        

	

        beta.curr.ci[i.var]<-
                (
                p.z^2*var.beta2+
                n/sigma.sq*(beta.2*dnorm(in.phi)*1/sigma.sq^.5)^2*var.betals+
                n/sigma.sq*(beta.2*dnorm(in.phi)*beta.ls*lambda.shrink[i.var]/(n-1))^2*var.sigma+
                n/var.betals*(beta.2*dnorm(in.phi)*ps.sigmasq/(n-1))^2*var.lambda
                )

	
	}
	beta.curr.ci<-abs(beta.curr.ci)^.5

	mean.temp<-beta.curr.mode
	set.temp<-beta.curr.mode==0
	if(sum(set.temp)>0)
        beta.curr.ci[set.temp]<-rnorm(sum(set.temp>0),mean.temp[set.temp],sd=beta.curr.ci[set.temp])
	set.temp<-beta.curr.mode>0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci[set.temp],
	lower=0,upper=Inf)
	set.temp<-beta.curr.mode<0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci[set.temp],
	lower=-Inf,upper=0)
	}
	#if(block){
	#A.big<-XprimeX+diag(1/Dtau)
	#A.big.inv<-ginv(A.big)
	#beta.var<-sigma.sq*A.big.inv
	#beta.vec<-A.big.inv%*%t(X.big)%*%(y-ran.effect)
	#	beta.curr<-mvrnorm(1,mu=beta.vec,Sigma=beta.var)
	#}


	
	if(k.raneff>0){
	re1<-updateREs(y, fits.main=as.vector(X.big%*%beta.curr), ran.effect1, ran.effect2, ran.effect3, id, id2, id3,	sigma.sq.b, sigma.sq.b2, sigma.sq.b3, 1,fix.eff=FALSE)
	sigma.sq.b<-re1$sigmas[1]
	sigma.sq.b2<-re1$sigmas[2]
	sigma.sq.b3<-re1$sigmas[3]
	ran.effect1<-re1$REs[,1]
	ran.effect2<-re1$REs[,2]
	ran.effect3<-re1$REs[,3]
	ran.effect<-rowSums(re1$REs)
	} 
	
	fits<-as.vector(X.big%*%beta.curr)+ran.effect
	##Update sigmas, lambda
	sigma.sq.shape<-(n-1)/2+p/2+k.raneff/2
	ps.sigma.sq.shape<-(n^.5-1)/2+p/2+k.raneff/2
	sigma.sq.scale<-sum((y-fits)^2)/2+ sum(beta.curr^2/Dtau)/2
	sigma.sq<-rinvgamma(1,shape=sigma.sq.shape,scale=sigma.sq.scale)
	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)
	if(type=="probit") sigma.sq<-1



		
	for(i in unique(c.cluster)) {
		lambda.all[c.cluster==i]<-rgamma(1, shape=sum(c.cluster==i)*n^.5+1, rate=sum(Dtau[c.cluster==i])/2+1.78)^.5
		}
	if(i.thin==1) p.adj<-n^.5
	if(one.constraint) lambda.all[1:length(lambda.all)]<-rgamma(1, shape=p*p.adj+1, rate=sum(Dtau)/2+1.78)^.5
	lambda.shrink<-lambda.all

	if(type=="probit"){
		y[y0==1]<-rtnorm(fits[y0==1],lower=0,upper=Inf)
		y[y0==0]<-rtnorm(fits[y0==0],lower=-Inf,upper=0)
	}
	
	if(type=="tobit"){
		fits.trunc<-X.big%*%beta.curr
		if(lower.trunc==TRUE) y[y==trunc]<-rtnorm(sum(y==trunc),lower=-Inf,upper=trunc,mean=fits.trunc[y==trunc],sd=sigma.sq^.5)
		if(lower.trunc==FALSE) y[y==trunc]<-rtnorm(sum(y==trunc),lower=trunc,upper=Inf,mean=fits.trunc[y==trunc],sd=sigma.sq^.5)
	}
	}#End thin loop
	
	if(i.gibbs>burnin){
		beta.mean[,i.gibbs-burnin]<-beta.curr
		beta.mode[,i.gibbs-burnin]<-beta.curr.mode
	sigma.sq.run[i.gibbs-burnin]<-sigma.sq
	lambda.run[,i.gibbs-burnin]<-lambda.all
	beta.ci[,i.gibbs-burnin]<-beta.curr.ci

	}
	
	if(i.gibbs%in%seq(1,burnin+gibbs,length=8)&save.temp){
		save(file="temp_sparsereg",beta.mode)
	}
	
	
	if(i.gibbs==1) cat("    0% of Burnin Completed:", i.gibbs,"/",burnin,"\n")	
	if(i.gibbs==floor(burnin/4)) cat("   25% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin/2)) cat("   50% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(3*burnin/4)) cat("   75% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs==floor(burnin)) cat("  100% of Burnin Completed:", i.gibbs,"/",burnin,"\n")
	if(i.gibbs == burnin+1) cat("Step 3 of 4: Burnin Completed; Saving Samples \n")
	if(i.gibbs==floor(burnin+1)) cat("    0% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+(gibbs)/4)) cat("   25% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+(gibbs)/2)) cat("   50% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+3*(gibbs)/4)) cat("   75% of Samples Gathered:", i.gibbs-burnin,"/",gibbs,"\n")
	if(i.gibbs==floor(burnin+gibbs)) cat("  100% of Samples Gathered:", i.gibbs-burnin,"/",burnin,"\n")

	##Gather things
	}#End gibbs loop
	cat("Step 4 of 4: Gathering Output\n")
	
		for(i.scale in 1:nrow(beta.mean)) {
			beta.mean[i.scale,]<-beta.mean[i.scale,]/scale.back[i.scale]*scale.y
			beta.mode[i.scale,]<-beta.mode[i.scale,]/scale.back[i.scale]*scale.y
			beta.ci[i.scale,]<-beta.ci[i.scale,]/scale.back[i.scale]*scale.y

			}
	matrix.names<-matrix(NA,nrow=3,ncol=ncol(X.big))
	rownames(matrix.names)<-c("X","treat1","treat2")

	for(i.scale in 1:ncol(X.big)) X.big[,i.scale]<-X.big[,i.scale]*scale.back[i.scale]


	row.names(beta.mode)<-gsub("_",": ",row.names(beta.mode))
	row.names(beta.mean)<-gsub("_",": ",row.names(beta.mean))

	output<-list("beta.mode"=as.mcmc(t(beta.mode)),"beta.mean"=as.mcmc(t(beta.mean)),"beta.ci"=as.mcmc(t(beta.ci)),
	"sigma.sq"=sigma.sq.run*scale.y^2,"X"=X.big,"y"=y,"lambda"=lambda.run,
	"varmat"=matrix.convert,"baseline"=baseline.vec,"modeltype"="onestage")
	class(output)<-"sparsereg"
	return(output)
}#End function loop
	
#sparsereg<-function(y, X,...) 	subgroup(y, X, treat=NULL, ...)	
#csts<-function(y, X, id=NULL, id2=NULL, id3=NULL, ...) 	subgroup(y, X,id=id,id2=id2,id3=id3,...)


