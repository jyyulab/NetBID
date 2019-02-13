#combRowEvid.2grps
#d: dataframe as below: (first column is the annotation)
#			geneSymbol  DMSO_1 DMSO_2 DMSO_3    DEX_1    DEX_2    DEX_3
#P0112072    PRKAR2B 0.50099 1.2108 1.0524 -0.34881 -0.13441 -0.87112
#P0112073    PRKAR2B 1.84579 2.0356 2.6025  1.62954  1.88281  1.29604
#comp as below:
# [1] 0 0 0 1 1 1
# Levels: 0 1
#d<-d[1,]
combRowEvid.2grps<-function(d,comp,
                            family=gaussian
                            ,method=c('MLE','Bayesian'),pooling=c('full','no','partial'), 
                            
                            n.iter=1000
                            , 
                            prior.V.scale=0.02
                            ,
                            prior.R.nu=1
                            , 
                            prior.G.nu=2
                            , 
                            nitt = 13000
                            , 
                            burnin =3000
                            ,
                            thin=10
                            ,
                            restand=TRUE
                            ,
                            logTransformed=TRUE
                            ,
                            log.base=2
                            ,
                            average.method=c('geometric','arithmetic')
                            ,
                            pseudoCount=0
){
  
  if(!all(comp %in% c(1,0))){
    stop('comp only takes 1 or 0 !!! \n')
  }
  
  if(missing(pooling))
    pooling<-c('full','no','partial')
  else
    pooling<-tolower(pooling)
  
  pooling<-match.arg(pooling,several.ok=T)
  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (!family$family %in% c('gaussian','binomial', 'poisson')) {
    print(family)
    stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
  }
  
  if(missing(method))
    method<-'Bayesian'
  if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
  if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'	
  method<-match.arg(method)
  
  
  #remove NAs
  nna<-apply(!is.na(d),2,all)
  d<-d[,nna]
  #comp<-comp[nna[-1]]
  
  #d<-d[1,]
  #cat(dim(d),'\n')
  #cat(d[1,1],'\n')
  d<-data.frame(t(d[,-1]))
  #cat(dim(d),'\n')
  
  
  
  
  dat<-data.frame(
    response=
      c(unlist(d[comp==levels(comp)[1],]),unlist(d[comp==levels(comp)[2],]))
    ,
    treatment=
      factor(c(rep(levels(comp)[1],sum(comp==levels(comp)[1])*ncol(d)),rep(levels(comp)[2],sum(comp==levels(comp)[2])*ncol(d))))
    ,
    probe=
      factor(c(rep(colnames(d),each=sum(comp==levels(comp)[1])),rep(colnames(d),each=sum(comp==levels(comp)[2]))))
  )
  
  
  #calculate FC
  FC.val<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='geometric',pseudoCount=pseudoCount)
  
  AveExp<-mean(dat$response)
  n.levels<-nlevels(dat$probe)
  
  rs<-c(FC=FC.val,AveExp=AveExp,n.levels=as.integer(n.levels))
  
  #cat(rs,'\n')
  
  if('arithmetic' %in% average.method){
    FC.ari<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='arithmetic',pseudoCount=0)	
    rs<-c(FC.ari_raw=FC.ari,rs)
  }
  
  
  #re-standarize the input data
  if(restand & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)
  
  #MLE approach to estimate parameters	
  if(method=='MLE'){
    if(family$family=='gaussian'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat)
        else
          M.no<-glm(response ~ treatment, data=dat)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      
      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){	
        if(n.levels>1){
          #consider sd==0 case
          if(sd(dat$response)==0)
          {
            dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
          }
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(	sum.tmp$devcomp$dims['n']-	sum.tmp$devcomp$dims['p']-	sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          #########################################
          #Another way to calculate pvalue
          #########################################
          #M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
          #anova(M.partial, M.null)
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))					
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)	
      }
    }else if(family$family=='binomial'){
      if('full'%in%pooling){
        M.full<-glm(treatment ~ response, data=dat, family=family)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      if('no'%in%pooling){	
        if(n.levels>1)
          M.no<-glm(treatment ~ response + probe, data=dat, family=family)
        else
          M.no<-glm(treatment ~ response, data=dat, family=family)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)	
      }
      
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            Dev.partial=as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(treatment ~ response, data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
      
      
    }else if(family$family=='poisson'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
        else
          M.no<-glm(response ~ treatment, data=dat,family='poisson')
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            #REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))	
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            #REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
      
    }
    else{
      stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
    }
  }
  
  #Bayesian approach
  else if(method=='Bayesian'){
    require(arm)
    require(MCMCglmm)
    
    if(family$family=='gaussian'){
      if('full'%in%pooling){
        #comlete pooling
        #				M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter)
        if(sd(dat$response)>0){	
          M.full<-bayesglm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.full)
          
          rs.full<-c(
            coef.full=sum.tmp$coef[2,1],
            se.full=sum.tmp$coef[2,2],
            t.full=sum.tmp$coef[2,3],
            pval.full=sum.tmp$coef[2,4],
            z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
            df.full=sum.tmp$df.residual-sum.tmp$df[1],
            AIC.full=sum.tmp$aic,
            BIC.full=AIC(M.full,k=log(nrow(dat))),
            Dev.full=sum.tmp$deviance
          )
        }else{
          rs.full<-c(
            coef.full=0,
            se.full=0,
            t.full=0,
            pval.full=1,
            z.full=0,
            df.full=NA,
            AIC.full=NA,
            BIC.full=NA,
            Dev.full=NA
          )
        }
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){				
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment + probe, data=dat)
        else
          #M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.no)
        
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
        
      }
      
      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))
        
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))
          
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          
          df.partial<-nrow(dat)-2
          
          #cat(n.levels,df.partial,'\n')
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        if(is.na(pval.partial))
          pval.partial<-sum.tmp$sol[2,5]
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
      
    }else if(family$family=='binomial'){
      
      if(family$link=='logit')
        prior.scale<-2.5
      else if(family$link=='probit')
        prior.scale<-2.5*1.6
      
      if('full'%in%pooling){
        #M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
        M.full<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #a bug in bayeglm for summary, fix: total obs - rank of the model
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){
        if(n.levels>1)
          #M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family,prior.scale = prior.scale)
        else
          #M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      if('partial'%in%pooling){	
        #categorial=logit, ordinal=probit
        if(family$link=='logit'){
          glmm.family<-'categorial'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels))					
          }
        }
        else if(family$link=='probit'){
          glmm.family<-'ordinal'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1))					
          }
        }
        else
          stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')
        
        #			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
        #			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
        #			scale<-var(dat$treatment)/var(dat$response)
        ################################# prior to be modified #################################
        #alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
        #########################################################################################
        sd.threshold<-prior.scale*4
        #coef.sign<- -sign(logFC)
        #				while(sd.threshold>prior.scale*2){				
        if(n.levels>1){
          M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        #coef.sign<-sign(rs.partial['coef.partial'])
        sd.threshold<-as.double(rs.partial['sd.partial'])
        #				}
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='poisson'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
        M.full<-bayesglm(response ~ treatment, data=dat,family='poisson')
        
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){				
        #no pooling
        if(n.levels>1)
          #M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment + probe, data=dat,family='poisson')
        else
          #M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))
        
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-(n.levels+1)*2
        }
        else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-2				
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else{
      stop('Only liner model and binomial family model are supported !!! \n')
    }
  }
  
  #taking care of extreme cases where se or sd is zero
  w1 <- grep('^t|^z',names(rs))
  w2 <- grep('^se|^sd',names(rs))
  w3 <- grep('^pval',names(rs))
  if(rs[w2]==0 & is.na(rs[w1[1]])==TRUE){
    rs[w1]<-0
    rs[w3]<-1
  }
  rs
}



#comparsion of mult-igroups
#only support linear Gausian model so far
combRowEvid.multgrps<-function(d,comp,family=gaussian,method=c('MLE','Bayesian'),pooling=c('full','no','partial'), n.iter=1000, prior.V.scale=0.02, prior.R.nu=1, prior.G.nu=2, nitt = 13000, burnin = 3000,thin=10,restand=TRUE,logTransformed=TRUE,log.base=2,pseudoCount=0){
  
  require(reshape)
  if(missing(pooling))
    pooling<-c('full','no','partial')
  else
    pooling<-tolower(pooling)
  
  pooling<-match.arg(pooling,several.ok=T)
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  #	if (!family$family %in% c('gaussian','binomial', 'poisson')) {
  if (!family$family %in% c('gaussian')) {		
    print(family)
    #stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
    stop("Only Gaussian model is supported!!! \n")		
  }
  
  if(missing(method))
    method<-'Bayesian'
  if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
  if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'	
  method<-match.arg(method)
  
  d<-reshape::melt(t(d[,-1]))
  
  dat<-data.frame(response=d$value,treatment=rep(comp,nlevels(d$X2)),probe=d$X2)
  
  AveExp<-mean(dat$response)
  n.levels<-nlevels(dat$probe)
  n.treatments<-nlevels(dat$treatment)
  rs<-c(AveExp=AveExp,n.levels=as.integer(n.levels))
  
  #re-standarize the input data
  if(restand & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)
  
  #MLE approach to estimate parameters	
  if(method=='MLE'){
    if(family$family=='gaussian'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ as.ordered(treatment), data=dat)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[-1,1],
          se.full=sum.tmp$coef[-1,2],
          t.full=sum.tmp$coef[-1,3],
          pval.full=sum.tmp$coef[-1,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }			
      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat)
        else
          M.no<-glm(response ~ treatment, data=dat)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2:n.treatments,1],
          se.no=sum.tmp$coef[2:n.treatments,2],
          t.no=sum.tmp$coef[2:n.treatments,3],
          pval.no=sum.tmp$coef[2:n.treatments,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      
      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){	
        if(n.levels>1){
          #consider sd==0 case
          if(sd(dat$response)==0)
          {
            dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
          }
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
          sum.tmp<-summary(M.partial)
          
          t.partial<-sum.tmp$coef[-1,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          #########################################
          #Another way to calculate pvalue
          #########################################
          #M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
          #anova(M.partial, M.null)
          rs.partial<-c(
            coef.partial=sum.tmp$coef[-1,1],
            se.partial=sum.tmp$coef[-1,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[-1,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))					
          rs.partial<-c(
            coef.partial=sum.tmp$coef[-1,1],
            se.partial=sum.tmp$coef[-1,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)	
      }
      
    }else if(family$family=='binomial'){
      if('full'%in%pooling){
        M.full<-glm(treatment ~ response, data=dat, family=family)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      if('no'%in%pooling){	
        if(n.levels>1)
          M.no<-glm(treatment ~ response + probe, data=dat, family=family)
        else
          M.no<-glm(treatment ~ response, data=dat, family=family)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)	
      }
      
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            Dev.partial=as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(treatment ~ response, data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
      
      
    }else if(family$family=='poisson'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
        else
          M.no<-glm(response ~ treatment, data=dat,family='poisson')
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            #REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))	
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            #REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
      
    }
    else{
      stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
    }
  }
  
  #Bayesian approach
  else if(method=='Bayesian'){
    require(arm)
    require(MCMCglmm)
    if(family$family=='gaussian'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter)
        M.full<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.full)
        
        rs.full<-c(
          coef.full=sum.tmp$coef[-1,1],
          se.full=sum.tmp$coef[-1,2],
          t.full=sum.tmp$coef[-1,3],
          pval.full=sum.tmp$coef[-1,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){				
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment + probe, data=dat)
        else
          #					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter)
          M.no<-bayesglm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2:n.treatments,1],
          se.no=sum.tmp$coef[2:n.treatments,2],
          t.no=sum.tmp$coef[2:n.treatments,3],
          pval.no=sum.tmp$coef[2:n.treatments,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(nlevels(dat$treatment))*prior.V.scale, nu = prior.G.nu )))
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*n.treatments
        }else{
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-n.treatments
        }
        sum.tmp<-summary(M.partial)
        sd.partial<-apply(M.partial$Sol[,-1],2,sd)
        t.partial<-sum.tmp$sol[-1,1]/sd.partial
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[-1,1],
          'l-95% CI.partial'=sum.tmp$sol[-1,2],
          'u-95% CI.partial'=sum.tmp$sol[-1,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[-1,5],
          #effective sample size
          n.effet.partial=sum.tmp$sol[-1,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd.partial,
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='binomial'){
      
      if(family$link=='logit')
        prior.scale<-2.5
      else if(family$link=='probit')
        prior.scale<-2.5*1.6
      
      if('full'%in%pooling){
        #M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
        M.full<-bayesglm(treatment ~ response, data=dat, family=family, prior.scale = prior.scale)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #a bug in bayeglm for summary, fix: total obs - rank of the model
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){
        if(n.levels>1)
          #M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family,prior.scale = prior.scale)
        else
          #					M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
          M.no<-bayesglm(treatment ~ response, data=dat, family=family,prior.scale = prior.scale)
        
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      if('partial'%in%pooling){	
        #categorial=logit, ordinal=probit
        if(family$link=='logit'){
          glmm.family<-'categorial'
          prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1))) 	
        }
        else if(family$link=='probit'){
          glmm.family<-'ordinal'
          prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1))) 	
        }
        else
          stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')
        
        #			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
        #			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
        #			scale<-var(dat$treatment)/var(dat$response)
        ################################# prior to be modified #################################
        #alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
        #########################################################################################
        sd.threshold<-prior.scale*4
        #coef.sign<- -sign(logFC)
        #				while(sd.threshold>prior.scale*2){				
        if(n.levels>1){
          M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        #coef.sign<-sign(rs.partial['coef.partial'])
        sd.threshold<-as.double(rs.partial['sd.partial'])
        #				}
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='poisson'){
      if('full'%in%pooling){
        #comlete pooling
        #M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
        M.full<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      
      if('no'%in%pooling){				
        #no pooling
        if(n.levels>1)
          #					M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment + probe, data=dat, family='poisson')
        else
          #					M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
          M.no<-bayesglm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      
      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu ))) 
        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-(n.levels+1)*2
        }
        else{
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-2				
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(			
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation 
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else{
      stop('Only liner model and binomial family model are supported !!! \n')
    }
  }
  rs
}
##
#fold change function
#positive: class1/class0
#negative: class0/class1
FC <-
  function(x,
           cl,
           logTransformed = TRUE,
           log.base = 2,
           average.method = c('geometric', 'arithmetic'),
           pseudoCount = 0) {
    x.class0 <- x[(cl == 0)] + pseudoCount
    x.class1 <- x[(cl == 1)] + pseudoCount
    if (missing(average.method))
      average.method <- 'geometric'
    if (logTransformed) {
      if (is.na(log.base) | log.base < 0)
        stop('You must specify log.bsae !\n')
      logFC <- mean(x.class1) - mean(x.class0)
      FC.val <- sign(logFC) * log.base ^ abs(logFC)
    } else{
      logFC <-
        ifelse(average.method == 'arithmetic',
               log(mean(x.class1)) - log(mean(x.class0)),
               mean(log(x.class1) - mean(log(x.class0))))
      FC.val <- sign(logFC) * exp(abs(logFC))
    }
    FC.val[FC.val == 0 | is.na(FC.val)] <- 1
    FC.val
  }