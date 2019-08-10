
############bid: Bayesian Inference of Drivers by differential analysis of case vs. ctrl
#' @export
bid<-function(eset,group.case,group.ctrl,signed=FALSE,
family=gaussian,
method='Bayesian',
n.iter=5000,
nitt=25000,
burnin=5000,
thin=1,
pooling='full',
logTransformed=TRUE,
restand=FALSE,
average.method='geometric'
,...){

  if(!group.case%in%pData(eset)$group)
    stop(paste('No group.case \'',group.case,'\' is defined in pData(eset)$group!\n',sep=''))


  if(!group.ctrl%in%pData(eset)$group)
    stop(paste('No group.ctrl \'',group.ctrl,'\' is defined in pData(eset)$group!\n',sep=''))


  if(all(fData(eset)[,1]==featureNames(eset))){
    names(fData(eset))[1]<-'id'
  }else{
    names(fData(eset))<-gsub('^id$','id.old',names(fData(eset)))
    fData(eset)<-data.frame(id=featureNames(eset),fData(eset))
  }


  rs<-fData(eset)

  eset.sel<-eset[,pData(eset)$group%in%c(group.case,group.ctrl)]

  comp <- factor(as.numeric(gsub(paste('^',group.ctrl,'$',sep=''),0,gsub(paste('^',group.case,'$',sep=''),1,(pData(eset.sel)$group)))));table(comp)

  d<-data.frame(id=featureNames(eset.sel),exprs(eset.sel),stringsAsFactors=FALSE)

  de<-plyr::ddply(d,'id','combRowEvid.2grps',comp=comp,family=family,method=method,n.iter=n.iter,nitt=nitt,burnin=burnin,thin=thin,pooling=pooling,logTransformed=logTransformed,restand=restand,average.method=average.method,...)

  names(de)<-gsub(paste('.',pooling,'$',sep=''),'',names(de))

  de$FDR.BH<-p.adjust(de$pval,'BH')

  de$log2FC<-sign(de$FC)*log2(abs(de$FC))


  if(signed) {
    de<-de[,c('id','z','pval','FDR.BH','log2FC')]
  }else{
    de$z<-abs(de$z)
    de<-de[,c('id','z','pval','FDR.BH')]
  }

  rs<-merge(rs,de,by='id',all.x=TRUE,sort=FALSE)

  dplyr::arrange(rs,pval)
}





#combine Pvalues by Fisher or Stouffer's method
#input:
#dat, a dataframe containing pvalue.cols
#sign.cols: column labels indicating sign of stat, if NULL, use absolute values to calculate overall stat for Stouffer's method
combinePvalues<-function(dat,pvalue.cols,sign.cols=NULL,FC.cols=NULL,logTransformed=FALSE,log.base=2,method=c('Stouffer','Fisher'),twosided=TRUE,signed=TRUE,byRow=FALSE){
  if(!is.data.frame(dat) & !is.matrix(dat))
    stop('dat must be in data.frame or matrix format!!! \n')

  dat.pvals<-as.matrix(dat[,pvalue.cols])

  if(missing(method))
    method<-'Stouffer'

  if(!is.null(sign.cols)){
    signs<-sign(as.matrix(dat[,sign.cols]))
    signs[signs==0]<-1
    if(ncol(signs) != ncol(dat.pvals)){
      stop('sign.cols must have the same length with pvalue.cols !!!')
    }
  }else{
    signs<-matrix(1,ncol=ncol(dat.pvals),nrow=nrow(dat.pvals))
  }
  if(byRow){
    rs0<-apply(dat.pvals*signs,1,combinePvalVector,method=method,twosided=twosided,signed=signed)
    rs<-t(rs0)
    colnames(rs)<-c('z','pval')
    nPvals<-unlist(apply(dat.pvals,1,f<-function(v){
      v<-as.vector(v);sum(!(is.na(v)|is.null(v)))
    }))
    if(!is.null(FC.cols)){
      FCs<-as.data.frame(dat[,FC.cols])
      if(logTransformed){
        logFC<-apply(FCs,1,mean,na.rm = TRUE)
      }else{
        FCs[FCs==0]<-1
        logFC<-apply(sign(FCs)*log(abs(FCs),base=log.base),1,mean,na.rm = TRUE)
        logFC<-sign(logFC)*(log.base^abs(logFC))
      }
      rs<-data.frame(nPvals=nPvals,FC=logFC,rs)
      names(rs)[2]<-ifelse(logTransformed,paste('log',log.base,'FC',sep=''),'FC')
    }else{
      rs<-data.frame(nPvals=nPvals,rs)
    }

    rs$FDR.BH<-p.adjust(rs$pval,method='BH')

  }else{
    rs0<-apply(dat.pvals*signs,2,combinePvalVector,method=method,twosided=twosided,signed=signed)
    rs<-as.vector(rs0)
    if(ncol(rs0)==1)
      names(rs)<-c('z','pval')
    else
      names(rs)<-paste(rep(c('z','pval'),ncol(rs0)),rep(colnames(rs0),each=2),sep='.')
    if(!is.null(FC.cols)){
      #		FC<-apply(as.data.frame(dat[,FC.cols]),2,mean,na.rm = TRUE)
      FCs<-as.data.frame(dat[,FC.cols])
      if(logTransformed){
        logFC<-apply(FCs,2,mean,na.rm = TRUE)
      }else{
        FCs[FCs==0]<-1
        logFC<-apply(sign(FCs)*log(abs(FCs),base=log.base),2,mean,na.rm = TRUE)
        logFC<-sign(logFC)*(log.base^abs(logFC))
      }
      if(length(logFC)==1)
        names(logFC)<-'log2FC'
      rs<-c(nPvals=nrow(dat.pvals),logFC,rs)
    }else{
      rs<-c(nPvals=nrow(dat.pvals),rs)
    }
  }
  rs
}


#Debugging
#pvals<-c(0.21,-.05)
#pvals<-c(runif(3,0,0.5),runif(3,-0.15,0)) ;pvals
#combinePvalVector(pvals,method='Fisher',signed=T,twosided = T)
#combinePvalVector(pvals,signed=T, twosided = T)
#combinePvalVector(c(1,1),method='Fisher')
#pvals<-c(1,1)
#combine one pvalues vector, the signed version
# For Fisher's method: if twosided, pvalues are transformed into single
combinePvalVector<-function(pvals,method=c('Stouffer','Fisher'),signed=TRUE,twosided=TRUE){

  #remove NA pvalues
  pvals<-pvals[!is.na(pvals) & !is.null(pvals)]

  if(sum(is.na(pvals))>=1){
    stat<-NA
    pval<-NA
  }else{
    if(twosided & (sum(pvals>1 | pvals< -1)>=1))
      stop('pvalues must between 0 and 1!\n')
    if(!twosided & (sum(pvals>0.5 | pvals< -0.5)>=1))
      stop('One-sided pvalues must between 0 and 0.5!\n')

    if(missing(method))
      method<-'Stouffer'

    if(!signed){
      pvals<-abs(pvals)
    }

    signs<-sign(pvals)
    signs[signs==0]<-1

    if(grepl('Fisher',method,ignore.case = TRUE)){
      if(twosided & signed){
        neg.pvals<-pos.pvals<-abs(pvals)/2
        pos.pvals[signs<0]<-1-pos.pvals[signs<0]
        neg.pvals[signs>0]<-1-neg.pvals[signs>0]
      }else{
        neg.pvals<-pos.pvals<-abs(pvals)
      }

      pvals<-c(1,-1)*c(pchisq(-2*sum(log(as.numeric(pos.pvals))),df=2*length(pvals),lower.tail = FALSE)/2,pchisq(-2*sum(log(as.numeric(neg.pvals))),df=2*length(pvals),lower.tail = FALSE)/2)

      pval<-min(abs(pvals))[1]
      #if two pvals are equal, pick up the first one
      stat <- sign(pvals[abs(pvals)==pval])[1]*qnorm(pval,lower.tail=F)[1]
      pval<-2*pval
    }
    else if(grepl('Stou',method,ignore.case = TRUE)){
      if(twosided){
        zs<-signs*qnorm(abs(pvals)/2,lower.tail=FALSE)
        stat<-sum(zs)/sqrt(length(zs))
        pval<-2*pnorm(abs(stat),lower.tail=FALSE)
      }
      else{
        zs<-signs*qnorm(abs(pvals),lower.tail=FALSE)
        stat<-sum(zs)/sqrt(length(zs))
        pval<-pnorm(abs(stat),lower.tail=FALSE)
      }
    }
    else{
      stop('Only \"Fisher\" or \"Stouffer\" method is supported!!!\n')
    }
  }
  return(c(stat=stat,pvalue=pval))
}








#d: dataframe as below: (first column is the annotation)
#			geneSymbol  DMSO_1 DMSO_2 DMSO_3    DEX_1    DEX_2    DEX_3
#P0112072    PRKAR2B 0.50099 1.2108 1.0524 -0.34881 -0.13441 -0.87112
#P0112073    PRKAR2B 1.84579 2.0356 2.6025  1.62954  1.88281  1.29604
#comp as below:
# [1] 0 0 0 1 1 1
# Levels: 0 1
#d<-d[1,]
#' @export
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
  comp<-comp[nna[-1]]

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
  if(is.na(rs[grep('^t|^z',names(rs))]) & rs[grep('^se|^sd',names(rs))]==0){
    rs[grep('^t|^z',names(rs))]<-0
    rs[grep('^pval',names(rs))]<-1
  }

  rs

}



