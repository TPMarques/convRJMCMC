CZ_ANOVA<-function(theta,chains,models,mcmciterations,nbatches=20,batchsize=max(mcmciterations)/nbatches,confidence=0.95,division="Batch"){

  mcmcdata<-data.frame(cbind(theta,chains,models,mcmciterations))

  colnames(mcmcdata)<-c("theta","chains","models","mcmciterations")

  if(division=="Batch"){

    mcmcdata$batch<-ceiling(mcmciterations/batchsize)
    mcmcdata_samples<-split(mcmcdata,mcmcdata$batch)

  }else{

    if(division=="Sequential"){

      mcmcdata_samples<-list()
      for(i in 1:nbatches){
        inf<-seq.int(0,max(mcmciterations)/2,by=(max(mcmciterations)/2)/nbatches)
        sup<-seq.int(max(mcmciterations)/2,max(mcmciterations),by=(max(mcmciterations)/2)/nbatches)
        mcmcdata_samples[[i]]<-with(mcmcdata,subset(mcmcdata,mcmciterations>inf[i]&mcmciterations<sup[i]))
      }

    }else{

      print("Undefined Splitting Method")

    }

  }




  result<-lapply(mcmcdata_samples,FUN=function(x){
    #Calculate ANOVAs
    anova1<-aov(theta~factor(chains),data=x)
    anova2<-aov(theta~factor(models),data=x)
    anova3<-aov(theta~factor(models)*factor(chains),data=x)

    #Get ANOVAs summary
    sm_aov1<-summary(anova1)
    sm_aov2<-summary(anova2)
    sm_aov3<-summary(anova3)

    #Calculate statistics of interest from ANOVAs summary
    V<-sum(sm_aov1[1][[1]][[2]])/sum(sm_aov1[1][[1]][[1]])
    Wc<-sm_aov1[1][[1]][[3]][[2]]
    Wm<-sm_aov2[1][[1]][[3]][[2]]
    Wcm<-sm_aov3[1][[1]][[3]][[4]]

    #Compute potential scale reduction factors and its confidence bounds
    PSRF1<-V/Wc
    ub_PSRF1<-PSRF1*qf(confidence,sm_aov1[1][[1]][[1]][[2]],sum(sm_aov1[1][[1]][[1]]))
    PSRF2<-Wm/Wcm
    ub_PSRF2<-PSRF2*qf(confidence,sm_aov3[1][[1]][[1]][[4]],sm_aov2[1][[1]][[1]][[2]])
    result<-list(sm_aov1,sm_aov2,sm_aov3,PSRF1,PSRF2,ub_PSRF1,ub_PSRF2,V,Wc,Wm,Wcm)
    names(result)<-c("Summary ANOVA 1","Summary ANOVA 2",
                     "Summary ANOVA 3","PSRF1","PSRF2",
                     "Upper Bound PSRF1","Upper Bound PSRF2",
                     "V","Wc","Wm","WcWm")
    return(result)
  })
  attr(result,"class")<-"CZ_ANOVA"
  names(result)<-"CZ_ANOVA"
  return(result)
}

CZ_MANOVA<-function(theta,chains,models,mcmciterations,nbatches=20,batchsize=max(mcmciterations)/nbatches,confidence=0.95,division="Batch"){

  mcmcdata<-data.frame(cbind(theta,chains,models,mcmciterations))

  thetanames<-c()

  for(i in 1:ncol(theta)){
    thetanames<-c(thetanames,paste("theta",i,sep=""))
  }

  colnames(mcmcdata)<-c(thetanames,"chains","models","mcmciterations")

  if(division=="Batch"){

    mcmcdata$batch<-ceiling(mcmciterations/batchsize)
    mcmcdata_samples<-split(mcmcdata,mcmcdata$batch)

  }else{

    if(division=="Sequential"){

      mcmcdata_samples<-list()
      for(i in 1:nbatches){
        inf<-seq.int(0,max(mcmciterations)/2,by=(max(mcmciterations)/2)/nbatches)
        sup<-seq.int(max(mcmciterations)/2,max(mcmciterations),by=(max(mcmciterations)/2)/nbatches)
        mcmcdata_samples[[i]]<-with(mcmcdata,subset(mcmcdata,mcmciterations>inf[i]&mcmciterations<sup[i]))
      }

    }else{

      print("Undefined Splitting Method")

    }

  }
  result<-lapply(mcmcdata_samples,FUN=function(x){
    MANOVA1<-manova(sapply(thetanames,FUN=function(y){with(x,get(y))})~factor(chains),data=x)
    MANOVA2<-manova(sapply(thetanames,FUN=function(y){with(x,get(y))})~factor(models),data=x)
    MANOVA3<-manova(sapply(thetanames,FUN=function(y){with(x,get(y))})~factor(models)*factor(chains),data=x)

    V<-(summary(MANOVA1)[2]$SS$'factor(chains)'+summary(MANOVA1)[2]$SS$Residuals)/sum(as.numeric(unlist(summary(MANOVA1)[4])[1:2]))

    Wc<-summary(MANOVA1)[2]$SS$Residuals/as.numeric(unlist(summary(MANOVA1)[4])[2])

    MPSRF1<-max(eigen(solve(Wc)%*%V)$values)

    Wm<-summary(MANOVA2)[2]$SS$Residuals/as.numeric(unlist(summary(MANOVA2)[4])[2])

    WmWc<-summary(MANOVA3)[2]$SS$Residuals/unlist(summary(MANOVA3)[4])[4]

    MPSRF2<-max(eigen(solve(WmWc)%*%Wm)$values)

    V_eig<-max(eigen(V)$values)
    Wc_eig<-max(eigen(Wc)$values)
    Wm_eig<-max(eigen(Wm)$values)
    WmWc_eig<-max(eigen(WmWc)$values)

    result<-list(V,Wc,Wm,WmWc,MPSRF1,MPSRF2,V_eig,Wc_eig,Wm_eig,WmWc_eig)
    return(result)
  })
  attr(result,"class")<-"CZ_MANOVA"
  return(result)
}

plot.CZ_ANOVA<-function(cza_obj){
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  PSRF1<-as.vector(unlist(sapply(cza_obj,"[",4)))
  PSRF2<-as.vector(unlist(sapply(cza_obj,"[",5)))
  ub_PSRF1<-as.vector(unlist(sapply(cza_obj,"[",6)))
  ub_PSRF2<-as.vector(unlist(sapply(cza_obj,"[",7)))
  V<-as.vector(unlist(sapply(cza_obj,"[",8)))
  Wc<-as.vector(unlist(sapply(cza_obj,"[",9)))
  Wm<-as.vector(unlist(sapply(cza_obj,"[",10)))
  Wcm<-as.vector(unlist(sapply(cza_obj,"[",11)))
  ylim<-c(0.99*min(PSRF1,PSRF2),1.01*max(ub_PSRF1,ub_PSRF2))
  xv<-seq(1,length(PSRF1))
  plot(xv,PSRF1,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,PSRF2,col=2)
  lines(xv,ub_PSRF1,lty=3)
  lines(xv,ub_PSRF2,col=2,lty=3)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,1,3,3),col = c(1,2,1,2),legend = c("PSRF1","PSRF2","Limite Superior PSRF1","Limite Superior PSRF2"),ncol=4)
  devAskNewPage(ask = TRUE)
  ylim<-c(0.99*min(V,Wc),1.01*max(V,Wc))
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(xv,V,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,Wc,lty=3)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,3),col = c(1,1),legend = c("V","Wc"),ncol=2)
  devAskNewPage(ask = TRUE)
  ylim<-c(0.99*min(Wm,Wcm),1.01*max(Wm,Wcm))
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(xv,Wm,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,Wcm,lty=3)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,3),col = c(1,1),legend = c("Wm","WmWc"),ncol=2)
}

plot.CZ_MANOVA<-function(cza_obj){
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  MPSRF1<-as.vector(unlist(sapply(cza_obj,"[",5)))
  MPSRF2<-as.vector(unlist(sapply(cza_obj,"[",6)))
  V<-as.vector(unlist(sapply(cza_obj,"[",7)))
  Wc<-as.vector(unlist(sapply(cza_obj,"[",8)))
  Wm<-as.vector(unlist(sapply(cza_obj,"[",9)))
  Wcm<-as.vector(unlist(sapply(cza_obj,"[",10)))
  ylim<-c(0.99*min(MPSRF1,MPSRF2),1.01*max(MPSRF1,MPSRF2))
  xv<-seq(1,length(MPSRF1))
  plot(xv,MPSRF1,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,MPSRF2,col=2)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,1),col = c(1,2),legend = c("MPSRF1","MPSRF2"),ncol=2)
  devAskNewPage(ask = TRUE)
  ylim<-c(0.99*min(V,Wc),1.01*max(V,Wc))
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(xv,V,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,Wc,lty=3)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,3),col = c(1,1),legend = c("Maior autovalor de V","Maior autovalor de Wc"),ncol=2)
  devAskNewPage(ask = TRUE)
  ylim<-c(0.99*min(Wm,Wcm),1.01*max(Wm,Wcm))
  layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(xv,Wm,type="l",ylim=ylim,xlab="Lote",ylab=paste("Crit",intToUtf8(0233),"rio",sep=""),main=NULL)
  lines(xv,Wcm,lty=3)
  par(mar=c(rep(0,4)))
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  legend("center",lty=c(1,3),col = c(1,1),legend = c("Maior autovalor de Wm","Maior autovalor de WmWc"),ncol=2)
}
