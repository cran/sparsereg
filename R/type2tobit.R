########################################
########################################
####The functions
########################################
########################################

#Takes type linear, tobit, probit

#scale.type: none, expand.treat, expand.X

type2tobit<-function(y1, y2=NULL, X1, X2=NULL, stage2, 
group=NULL,  gibbs=200, burnin=200, 
thin=10, one.constraint=TRUE, id=NULL, id2=NULL, id3=NULL, 
id_2=NULL, id2_2=NULL, id3_2=NULL){


	
cat("Step 1 of 4: Formatting Data","\n")
	###################################
	###################################
	#### Declare variables
	###################################
	###################################
	cluster<-group
	set.group<-(length(group)==0)
	if(set.group){
		cluster<-group<-rep(1,ncol(X1))
		}
	if(length(group)!=ncol(X1)|(length(group)!=ncol(X2))) stop("group should have one entry for every element of the full X matrix")

	#Drop zero observations from second stage
	X2<-X2[y2>0,]
	if(length(id_2)>0) id_2<-id_2[y2>0]
	if(length(id2_2)>0) id2_2<-id2_2[y2>0]
	if(length(id3_2)>0) id3_2<-id3_2[y2>0]
	y2<-y2[y2>0]

	n<-length(y1)
	n_2<-length(y2)
	sa<-stage2

	if(length(cluster)==0) cluster<-rep(1,ncol(X2))
	c.cluster<-cluster<-group
	
	if(length(id)>0) id<-as.factor(as.numeric(as.factor(id)))
	if(length(id2)>0) id2<-as.factor(as.numeric(as.factor(id2)))
	if(length(id3)>0) id3<-as.factor(as.numeric(as.factor(id3)))

	if(length(id_2)>0)  id_2<-as.factor(as.numeric(as.factor(id_2)))
	if(length(id2_2)>0)  id2_2<-as.factor(as.numeric(as.factor(id2_2)))
	if(length(id3_2)>0)  id3_2<-as.factor(as.numeric(as.factor(id3_2)))

	k.effect<-length((unique(id))) + length(unique(id2)) + length(unique(id3))
	k.effect_2<-length((unique(id_2))) + length(unique(id2_2)) + length(unique(id3_2))



	###################################
	###################################
	#### Process data
	###################################
	###################################

	if(length(colnames(X1))==0) colnames(X1)<-paste("X1",1:ncol(X1),sep="_")
	if(length(colnames(X2))==0) colnames(X2)<-paste("X2",1:ncol(X2),sep="_")
	colnames.orig<-colnames(X1)
	colnames.orig_2<-colnames(X2)
		
	X1<-X1[,apply(X2,2,sd)>0]
	X2<-X2[,apply(X2,2,sd)>0]
	X.big<-X1	
	X.big<-as.matrix(X.big)
	X.big<-apply(X.big,2,as.numeric)
	scales.X<-apply(X.big,2,FUN=function(x) sum(x^2)^.5)
		X.big<-apply(X.big,2,FUN=function(x) x-mean(x) )
		keeps<-apply(X.big,2,sd)>0&(apply(X.big,2,FUN=function(x) sum(is.na(x)))==0)
		X.big<-X.big[,keeps]
		c.cluster<-c.cluster[keeps]

		cor1<-cor(X.big)
		diag(cor1)<-0
		cors.run<-NULL
	for(i in 1:(ncol(X.big)-1)) for(j in (i+1):ncol(X.big)){
		if(cor1[i,j]^2>0.9999) cors.run<-rbind(cors.run,c(i,j))
	}
		if(length(cors.run)>0){
		drops<-as.vector(cors.run[,2])
		X.big<-X.big[,-drops]
		c.cluster<-c.cluster[-drops]
		}

	X0<-X.big
	scale.back<-apply(X.big,2,sd)
	for(i.scale in 1:ncol(X.big)) X.big[,i.scale]<-X.big[,i.scale]/scale.back[i.scale]
	

	##Stage 1 dv
	z<-rep(dnorm(1)/pnorm(1),n)
	z[y1<=0] <- -z[y1<=0] 

	##Stage 2 dv and X matrix
	y_2<-y2
	X.big_2<-X2
	X.big_2<-apply(X.big_2,2,FUN=function(x) x-mean(x))

	rm(X0);rm(X1);rm(X2)
	gc()
	
	scale.back_2<-apply(X.big_2,2,sd)
	for(i.scale in 1:ncol(X.big_2)) X.big_2[,i.scale]<-X.big_2[,i.scale]/scale.back_2[i.scale]

	X.big_2<-X.big_2[,scale.back_2!=0]	
	scale.y_2<-sd(y_2)
	y_2<-(y_2-mean(y_2))/scale.y_2
	c.cluster_2<-c.cluster

	n.cluster<-length(unique(c.cluster))	
	n.cluster<-max(c.cluster)

	n.cluster_2<-length(unique(c.cluster_2))	
	n.cluster_2<-max(c.cluster_2)

	###################################
	###################################
	#### Format Containers and Initialize
	###################################
	###################################
	
	p<-ncol(X.big)
	p_2<-ncol(X.big_2)

	c.cluster<-c.cluster[1:p]
	c.cluster[!is.finite(c.cluster)]<-1
	c.cluster_2<-c.cluster_2[1:p_2]
	c.cluster_2[!is.finite(c.cluster_2)]<-1
	
	beta.mean<-matrix(NA,nrow=p, ncol=gibbs)
	beta.mode<-matrix(NA,nrow=p, ncol=gibbs)
	beta.ci<-matrix(NA,nrow=p,ncol=gibbs)

	beta.mean_2<-matrix(NA,nrow=p_2, ncol=gibbs)
	beta.mode_2<-matrix(NA,nrow=p_2, ncol=gibbs)
	beta.ci_2<-matrix(NA,nrow=p_2,ncol=gibbs)

	helpman.run<-matrix(NA,nrow=2, ncol=gibbs)
	lambda.run<-matrix(NA,nrow=p, ncol=gibbs)
	lambda.run_2<-matrix(NA,nrow=p_2, ncol=gibbs)

	sigma.sq.run<-rep(NA,gibbs)
	rownames(beta.mean)<-rownames(beta.mode)<-colnames(X.big)
	rownames(beta.mean_2)<-rownames(beta.mode_2)<-colnames(X.big_2)
	
	beta.curr.mode<-beta.curr<-rnorm(p,sd=.1)
	beta.curr.mode_2<-beta.curr_2<-rnorm(p_2,sd=.1)


	ps.sigmasq_2<-sigma.sq_2<-sigma.sq<-1#mean((y-X.big%*%beta.curr)^2)

	lambda.all<-rep(0,p)
	for(i.cluster in unique(c.cluster)) {lambda.all[c.cluster==i.cluster]<-(2*sum(c.cluster==i.cluster)/sum(abs(beta.curr)[c.cluster==i.cluster]))^.5
		}
	lambda.all[grep("intercept",colnames(X.big))]<-1e-4
	lambda.all_2<-lambda.all[1:length(beta.curr_2)]
	lambda.all_2[!is.finite(lambda.all_2)]<-max(lambda.all_2,na.rm=TRUE)

	k.vec<-rep(NA,length(c.cluster))
	for(i in 1:n.cluster) k.vec[c.cluster==i]<-sum(c.cluster==i)

	###################################
	###################################
	#### Initialize Random Effects
	###################################
	###################################
	
	ran.effect<-ran.effect1<-ran.effect2<-ran.effect3<-rep(0,n)
	sigma.sq.b<-sigma.sq.b2<-sigma.sq.b3<-1

	ran.effect_2<-ran.effect1_2<-ran.effect2_2<-ran.effect3_2<-rep(0,n_2)
	sigma.sq.b_2<-sigma.sq.b2_2<-sigma.sq.b3_2<-1
	select.adj_2<-rep(0,n_2)
	
	if(k.effect>0){
	re1<-updateREs(z, fits.main=as.vector(X.big%*%beta.curr), ran.effect1, ran.effect2, ran.effect3, id, id2, id3,sigma.sq.b, sigma.sq.b2,sigma.sq.b3,sigma.sq=1,fix.eff=FALSE)
	ran.effect<-rowSums(re1$REs)
	}
	
	if(k.effect_2>0){
	re2<-updateREs(y_2, fits.main=as.vector(X.big_2%*%beta.curr_2+select.adj_2),ran.effect1_2,ran.effect2_2,ran.effect3_2,id_2,id2_2,id3_2,	sigma.sq.b_2,sigma.sq.b2_2,sigma.sq.b3_2,sigma.sq_2,fix.eff=FALSE)
	ran.effect_2<-rowSums(re2$REs)
	}
	##Delta term from Helpman et al.
	helpman.part<-delta<-beta.mills<-0


	##Begin gibbs loop
	for(i.gibbs in 1:(burnin+gibbs)){	
	for(i.thin in 1:thin){#Start thin loop

	###################################
	###################################
	#### Stage 1: Probit
	###################################
	###################################

	##Update tau^2; need these for betas
	sigma.sq<-1
	beta.curr<-as.vector(beta.curr)
	muprime.all<-abs(lambda.all*sqrt(sigma.sq))/abs(beta.curr)
	lambda.prime<-lambda.all^2
		invTau2<-sapply(1:length(muprime.all), FUN=function(i) rinv.gaussian(1, muprime.all[i], (lambda.prime[i] ) ) )	
	Dtau<-abs(1/invTau2)
	mean.z<-mean(z)
	
	up1<-updatebeta_cpp(X0=X.big, y0=as.matrix(z-ran.effect-mean.z), betacurr0=as.matrix(beta.curr), betamode0=as.matrix(beta.curr.mode), lambdavec0=as.matrix(lambda.all), dtau0=as.matrix(Dtau), sigmasq0=as.matrix(1),ps_sigmasq0=as.matrix(n^.5),lambdashrink0=as.matrix(lambda.all),k0=as.matrix(1)  )

	beta.curr<-up1$beta.mean
	beta.curr.mode<-up1$beta.mode

	beta.curr.ci<-up1$beta.ci
	if(i.thin==thin){

#	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)

	for(i.var in 1:p){
        beta.ls<-up1$beta.ols[i.var]
        lambda.use<-lambda.all[i.var]
	in.phi<-abs(n^.5*abs(beta.ls/sigma.sq^.5)-lambda.all[i.var]*(n^.5)^.5/(sigma.sq^.5*(n-1)) ) 
	p.z<-pnorm(in.phi)
        beta.2<-beta.curr.mode[i.var]
	var.beta2<-sigma.sq/(n-1)
        var.betals<-beta.curr.ci[i.var]
	var.lambda<-1/(4*lambda.use^2)*(p*n^.5+1)/(sum(Dtau)/2+1.78)^2
        #beta.sigma<-sigma.sq.scale#(beta.curr.ci[i.var]*(n-1)+sum(beta.curr^2/Dtau))/2
        #alpha.sigma<-ps.sigma.sq.shape#(n^.5-1)/2
        #Variance of pseudo sigma
        var.sigma<-0#1/(4*ps.sigmasq)*beta.sigma/((alpha.sigma-1)^2*(alpha.sigma-2))
        

	

        beta.curr.ci[i.var]<-
                (
                p.z^2*var.beta2+
                n/sigma.sq*(beta.2*dnorm(in.phi)*1/sigma.sq^.5)^2*var.betals+
                n/sigma.sq*(beta.2*dnorm(in.phi)*beta.ls*lambda.all[i.var]/(n-1))^2*var.sigma+
                n/var.betals*(beta.2*dnorm(in.phi)*n^.5/(n-1))^2*var.lambda
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

	
		
	fits.main<-as.vector(X.big%*%beta.curr)
	if(k.effect>0){
	re1<-updateREs(z, fits.main=fits.main, ran.effect1, ran.effect2, ran.effect3, id, id2, id3,	sigma.sq.b, sigma.sq.b2, sigma.sq.b3, 1,fix.eff=FALSE)
	sigma.sq.b<-re1$sigmas[1]
	sigma.sq.b2<-re1$sigmas[2]
	sigma.sq.b3<-re1$sigmas[3]
	sigma.sq<-1
	ran.effect1<-re1$REs[,1]
	ran.effect2<-re1$REs[,2]
	ran.effect3<-re1$REs[,3]
	ran.effect<-rowSums(re1$REs)
	}

	theta.probit<-X.big%*%beta.curr+ran.effect+mean(z)
	z[y1>0]<-rtnorm(sum(y1>0),mean=theta.probit[y1>0],sd=1,lower=0,upper=Inf)
	z[y1==0]<-rtnorm(sum(y1==0),mean=theta.probit[y1==0],sd=1,lower=-Inf,upper=0)
	z[abs(z)>100]<-100*sign(z[abs(z)>100])

	for(i in 1:n.cluster) lambda.all[c.cluster==i]<-rgamma(1, shape=sum(c.cluster==i)*n^.5+1, rate=sum(Dtau[c.cluster==i])/2+1.78)^.5
	

	##Generate Inverse Mills ratio
	inv.mills<-rep(0,n)
	inv.mills<-dnorm(z,log=TRUE)-pnorm(z,log.p=TRUE)
	inv.mills<-exp(inv.mills)
	inv.mills[is.na(inv.mills)]<- min(inv.mills,na.rm=TRUE)/2
	inv.mills[is.infinite(inv.mills)]<-max(inv.mills[is.finite(inv.mills)],na.rm=TRUE)
	range(z)
	if(sum(sa)>0) inv.mills<-c(rep(0,sum(sa)),inv.mills[y1>0])
	if(sum(sa)==0) inv.mills<-inv.mills[y1>0]
	inv.mills<-inv.mills-mean(inv.mills)

	###################################
	###################################
	#### Stage 2 Updates
	###################################
	###################################
	
	
	beta.curr_2<-as.vector(beta.curr_2)
	muprime.all_2<-abs(lambda.all_2*sqrt(sigma.sq_2))/abs(beta.curr_2)
	lambda.prime_2<-lambda.all_2^2
		invTau2_2<-sapply(1:length(muprime.all_2), FUN=function(i) rinv.gaussian(1, muprime.all_2[i], (lambda.prime_2[i] ) ) )	
	Dtau_2<-abs(1/invTau2_2)

	up1_2<-updatebeta_cpp(X0=X.big_2, y0=as.matrix(y_2-ran.effect_2-helpman.part), betacurr0=as.matrix(beta.curr_2), betamode0=as.matrix(beta.curr.mode_2), lambdavec0=as.matrix(lambda.all_2), dtau0=as.matrix(Dtau_2), sigmasq0=as.matrix(sigma.sq_2), ps_sigmasq0=as.matrix(ps.sigmasq_2), lambdashrink0=as.matrix(lambda.all_2), k0=as.matrix(1)  )

	beta.curr_2<-up1_2$beta.mean
	beta.curr.mode_2<-up1_2$beta.mode
	beta.curr_2[!is.finite(beta.curr_2)]<-0
	beta.curr.mode_2[!is.finite(beta.curr.mode_2)]<-0

	beta.curr.ci_2<-up1_2$beta.ci
	if(i.thin==thin){

#	ps.sigmasq<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)

	for(i.var in 1:p){
        beta.ls<-up1_2$beta.ols[i.var]
        lambda.use_2<-lambda.all_2[i.var]
	in.phi<-abs(n^.5*abs(beta.ls/sigma.sq_2^.5)-lambda.all_2[i.var]*ps.sigmasq_2^.5/(sigma.sq_2^.5*(n-1)) ) 
	p.z<-pnorm(in.phi)
        beta.2<-beta.curr.mode_2[i.var]
	var.beta2<-sigma.sq_2/(n-1)
        var.betals<-beta.curr.ci_2[i.var]
	var.lambda_2<-1/(4*lambda.use_2^2)*(p*n^.5+1)/(sum(Dtau_2)/2+1.78)^2
        beta.sigma<-sigma.sq.scale#(beta.curr.ci[i.var]*(n-1)+sum(beta.curr^2/Dtau))/2
        alpha.sigma<-ps.sigma.sq.shape#(n^.5-1)/2
        #Variance of pseudo sigma
        var.sigma<-1/(4*ps.sigmasq_2)*beta.sigma/((alpha.sigma-1)^2*(alpha.sigma-2))
        

	

        beta.curr.ci_2[i.var]<-
                (
                p.z^2*var.beta2+
                n/sigma.sq*(beta.2*dnorm(in.phi)*1/sigma.sq_2^.5)^2*var.betals+
                n/sigma.sq*(beta.2*dnorm(in.phi)*beta.ls*lambda.all_2[i.var]/(n-1))^2*var.sigma+
                n/var.betals*(beta.2*dnorm(in.phi)*ps.sigmasq_2/(n-1))^2*var.lambda_2
                )

	
	}
	beta.curr.ci_2<-abs(beta.curr.ci_2)^.5

	mean.temp<-beta.curr.mode_2
	set.temp<-beta.curr.mode_2==0
	if(sum(set.temp)>0)
        beta.curr.ci[set.temp]<-rnorm(sum(set.temp>0),mean.temp[set.temp],sd=beta.curr.ci_2[set.temp])
	set.temp<-beta.curr.mode_2>0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci_2[set.temp],
	lower=0,upper=Inf)
	set.temp<-beta.curr.mode_2<0
	if(sum(set.temp)>0)
	beta.curr.ci[set.temp]<-rtnorm(sum(set.temp),mean=mean.temp[set.temp],sd=beta.curr.ci_2[set.temp],
	lower=-Inf,upper=0)
	}


	fits.main_2<-as.vector(X.big_2%*%beta.curr_2)
	if(k.effect_2>0){
	re2<-updateREs(y_2, fits.main=fits.main_2-helpman.part, ran.effect1_2, ran.effect2_2, ran.effect3_2, id_2, id2_2, id3_2,	sigma.sq.b_2, sigma.sq.b2_2, sigma.sq.b3_2, sigma.sq_2,fix.eff=FALSE)
	sigma.sq.b_2<-re2$sigmas[1]
	sigma.sq.b2_2<-re2$sigmas[2]
	sigma.sq.b3_2<-re2$sigmas[3]
	ran.effect1_2<-re2$REs[,1]
	ran.effect2_2<-re2$REs[,2]
	ran.effect3_2<-re2$REs[,3]
	ran.effect_2<-rowSums(re2$REs)
	}
	k.raneff_2<- length(unique(id_2))+length(unique(id2_2))+length(unique(id3_2))

	fits_2<-as.vector(X.big_2%*%beta.curr_2)+ran.effect_2+helpman.part
	sigma.sq.shape<-(n_2-1)/2+p_2/2+k.raneff_2/2
	ps.sigma.sq.shape<-(n_2^.5-1)/2+p/2+k.raneff_2/2
	sigma.sq.scale<-sum((y_2-fits_2)^2)/2+ sum(beta.curr_2^2/Dtau_2)/2

	sigma.sq_2<-rinvgamma(1,shape=sigma.sq.shape,scale=sigma.sq.scale)
	ps.sigmasq_2<-rinvgamma(1,shape=ps.sigma.sq.shape,scale=sigma.sq.scale)
		
	for(i in unique(c.cluster_2)) {
		lambda.all_2[c.cluster_2==i]<-rgamma(1, shape=sum(c.cluster_2==i)*n_2^.5+1, rate=sum(Dtau_2[c.cluster_2==i])/2+1.78)^.5
		}

	lambda.shrink_2<-lambda.all_2

	##Update helpman part
	lm1<-lm(y_2-fits_2+helpman.part~inv.mills)
	s1<-summary(lm1)
	beta.mills<-rnorm(1,mean=lm1$coef[2],sd=s1$coef[2,2])
	helpman.part<-inv.mills*beta.mills
	helpman.part<-helpman.part-mean(helpman.part)

	}#End thin loop

	###################################
	###################################
	#### Gather results
	###################################
	###################################

	if(i.gibbs>burnin){
		beta.mean[,i.gibbs-burnin]<-beta.curr
		beta.mode[,i.gibbs-burnin]<-beta.curr.mode
		beta.ci[,i.gibbs-burnin]<-beta.curr.ci
		beta.mean_2[,i.gibbs-burnin]<-beta.curr_2
		beta.mode_2[,i.gibbs-burnin]<-beta.curr.mode_2
		beta.ci_2[,i.gibbs-burnin]<-beta.curr.ci

		helpman.run[,i.gibbs-burnin]<-c(delta,beta.mills)

	sigma.sq.run[i.gibbs-burnin]<-sigma.sq_2		

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



		scale.back_2<-scale.back_2[scale.back_2!=0]
		same.scale<-FALSE
		if(same.scale==FALSE){
		for(i.scale in 1:nrow(beta.mean)) {
			beta.mean[i.scale,]<-beta.mean[i.scale,]/scale.back[i.scale]
			beta.mode[i.scale,]<-beta.mode[i.scale,]/scale.back[i.scale]
			beta.ci[i.scale,]<-beta.ci[i.scale,]/scale.back[i.scale]

			}
			for(i.scale in 1:nrow(beta.mean_2)) {
			beta.mean_2[i.scale,]<-beta.mean_2[i.scale,]/scale.back_2[i.scale]*scale.y_2
			beta.mode_2[i.scale,]<-beta.mode_2[i.scale,]/scale.back_2[i.scale]*scale.y_2
			beta.ci_2[i.scale,]<-beta.ci_2[i.scale,]/scale.back_2[i.scale]*scale.y_2

			}
		} else{
			beta.mean<-beta.mean/(nrow(X.big)-1)^.5
			beta.mode<-beta.mode/(nrow(X.big)-1)^.5
			}

	beta.mean.out<-beta.mode.out<-beta.mean_2.out<-beta.mode_2.out<-matrix(0,nrow=length(colnames.orig),ncol=ncol(beta.mean))
		rownames(beta.mean.out)<-rownames(beta.mode.out)<-rownames(beta.ci)<-colnames.orig
		rownames(beta.mean_2.out)<-rownames(beta.mode_2.out)<-rownames(beta.ci_2)<-colnames.orig_2
		
	for(i.name in 1:nrow(beta.mean)){
		beta.mean.out[rownames(beta.mean.out)==rownames(beta.mean)[i.name],]<-beta.mean[i.name,]
		}
	for(i.name in 1:nrow(beta.mode)){
		beta.mode.out[rownames(beta.mode.out)==rownames(beta.mode)[i.name],]<-beta.mode[i.name,]
		}	
	for(i.name in 1:nrow(beta.ci)){
		beta.ci[rownames(beta.ci)==rownames(beta.ci)[i.name],]<-beta.ci[i.name,]
		}	
	for(i.name in 1:nrow(beta.mean_2)){
		beta.mean_2.out[rownames(beta.mean_2.out)==rownames(beta.mean_2)[i.name],]<-beta.mean_2[i.name,]
		}		
	for(i.name in 1:nrow(beta.mode_2)){
		beta.mode_2.out[rownames(beta.mode_2.out)==rownames(beta.mode_2)[i.name],]<-beta.mode_2[i.name,]
		}
	for(i.name in 1:nrow(beta.ci_2)){
		beta.ci_2[rownames(beta.ci_2)==rownames(beta.ci_2)[i.name],]<-beta.ci_2[i.name,]
		}	
		

	output<-list("beta.mode"=as.mcmc(t(beta.mode.out)),"beta.mean"=as.mcmc(t(beta.mean.out)),
"beta.ci"=as.mcmc(t(beta.ci)),
	"beta.mode.2"=as.mcmc(t(beta.mode_2.out)),"beta.mean.2"=as.mcmc(t(beta.mean_2.out)), 
	"beta.ci_2"=as.mcmc(t(beta.ci_2)), "heckman.run"=helpman.run,
	"sigma.sq"=sigma.sq.run*scale.y_2^2,"X"=X.big,"modeltype"="twostage")
	class(output)<-"sparsereg"
	return(output)
}
	
