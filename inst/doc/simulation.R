## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,eval=FALSE)
#setwd("~/Desktop/cornet")
#devtools::install_github("rauschenberger/cornet")

## ----abstract-----------------------------------------------------------------
#  grDevices::pdf("manuscript/figure_idea.pdf",width=5,height=2.5)
#  
#  box <- function(x,y,width=0.22,height=0.2,labels="",cex=1,col="black",...){
#    xs <- x + 0.5*c(-1,-1,1,1)*width
#    ys <- y + 0.5*c(-1,1,1,-1)*height
#    graphics::polygon(x=xs,y=ys,border=col,lwd=2,...)
#    graphics::text(x=x,y=y,labels=labels,col=col,cex=cex)
#  }
#  
#  graphics::par(mar=c(0,0,0,0))
#  graphics::plot.new()
#  graphics::plot.window(xlim=c(0,1),ylim=c(0,1))
#  
#  v <- h <- 0.1
#  
#  box(x=0+h,y=0.5,labels="outcomes,\nfeatures")
#  box(x=0.5,y=1-v,labels="initial binary\nclassification",col="red")
#  box(x=0.5,y=0+v,labels="numerical\nprediction",col="blue")
#  box(x=1-h,y=0.5,labels="final binary\nclassification",col="red")
#  
#  d <- 0.02
#  graphics::arrows(x0=0.2+d,y0=0.5+c(-d,d),x1=0.4-d,y1=c(v,1-v),lwd=2,col=c("blue","red"))
#  graphics::arrows(x0=0.6+d,y0=c(v,1-v),x1=0.8-d,y1=0.5+c(-d,d),lwd=2,col=c("blue","red"))
#  
#  graphics::text(x=0.4,y=0.55,labels="binary outcome:\nlogistic regression",col="red",cex=0.7,pos=3)
#  graphics::text(x=0.4,y=0.45,labels="numerical outcome:\nlinear regression",col="blue",cex=0.7,pos=1)
#  graphics::text(x=0.63,y=0.5,labels="combine\npredicted\nprobabilities",col="darkgrey",cex=0.7)
#  graphics::text(x=0.8,y=0.3,labels="transform\npredicted values to\npredicted probabilities",col="darkgrey",cex=0.7,pos=1)
#  
#  grDevices::dev.off()

## ----examples-----------------------------------------------------------------
#  loss <- list()
#  for(i in seq_len(4)){
#    loss[[i]] <- list()
#    cat("mode:",i,"\n")
#    for(j in seq_len(100)){
#      set.seed(j)
#      cat("iteration:",j,"\n")
#      n0 <- 100; n1 <- 10000; p <- 500
#      n <- n0 + n1
#      X <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
#      beta <- stats::rbinom(n=p,size=1,prob=0.05)*stats::rnorm(n=p)
#      eta <- X %*% beta
#      epsilon <- stats::rnorm(n=n)
#      if(i==1){
#        y <- eta + epsilon
#      } else if(i==2){
#        y <- ifelse(eta<0,-2,+2)+epsilon
#        table(y>=0,eta>=0)
#      } else if(i==3){
#        y <- ifelse(eta<0,-sqrt(abs(eta+epsilon)),(eta+epsilon)^2)
#      } else if(i==4){
#        y <- eta + epsilon + stats::rbinom(n=n,size=1,prob=0.05)*(2*stats::rbinom(n=n,size=1,prob=0.5)-1)*1.5*max(abs(eta))
#      }
#      foldid <- rep(c(0,1),times=c(n0,n1))
#      loss[[i]][[j]] <- cornet::cv.cornet(y=y,cutoff=0,X=X,foldid.ext=foldid)
#    }
#  }
#  
#  save(loss,file="results/simulation.RData")
#  writeLines(text=capture.output(utils::sessionInfo(),cat("\n"),
#          sessioninfo::session_info()),con="results/info_sim.txt")

## ----examples_figure----------------------------------------------------------
#  load("results/simulation.RData")
#  
#  grDevices::pdf("manuscript/figure_EXA.pdf",width=5,height=5)
#  
#  graphics::par(mfrow=c(2,2),mar=c(2,2,1,1))
#  pos <- c(binomial=1,combined=2,gaussian=3)
#  col <- c(binomial="red",combined="grey",gaussian="blue")
#  cex <- 0.7
#  names <- c("binomial","combined","gaussian")
#  for(i in seq_len(4)){
#    frame <- as.data.frame(t(sapply(loss[[i]],function(x) x$deviance)))
#    graphics::boxplot(x=frame[,names],at=pos[names],col=col[names],cex.axis=cex,main=paste0("example ",i),cex.main=cex,axes=FALSE)
#    graphics::box()
#    graphics::axis(side=1,at=pos[names],labels=names,cex.axis=cex,tick=FALSE,line=-1)
#    graphics::axis(side=2,cex.axis=cex)
#    for(j in c("binomial","combined","gaussian")){
#      mean <- mean(frame[[j]])
#      graphics::points(x=pos[j],y=mean,pch=21,col="white",bg="black")
#      if(j=="combined"){next}
#      pvalue <- stats::wilcox.test(x=frame$combined,y=frame[[j]],alternative="less")$p.value
#      signif <- ifelse(pvalue<=0.05/8,"*","")
#      graphics::text(x=mean(c(pos["combined"],pos[[j]])),y=min(frame),labels=paste0("p=",format(pvalue,digits=2,scientific=TRUE),signif),pos=3,cex=0.7)
#    }
#  }
#  
#  grDevices::dev.off()

## ----analysis,eval=FALSE------------------------------------------------------
#  iter <- 1000
#  set.seed(1)
#  frame <- data.frame(cor=runif(n=iter,min=0,max=0.9),
#                      n=round(runif(n=iter,min=100,max=200))+10000,
#                      prob=runif(n=iter,min=0.01,max=0.1),
#                      sd=runif(n=iter,min=1,max=2),
#                      exp=runif(n=iter,min=0.1,max=2),
#                      frac=runif(n=iter,min=0.5,max=0.9))
#  
#  ridge <- lasso <- list()
#  pb <- utils::txtProgressBar(min=0,max=nrow(frame),width=20,style=3)
#  for(i in seq_len(nrow(frame))){
#      utils::setTxtProgressBar(pb=pb,value=i)
#      set.seed(i)
#      data <- do.call(what=cornet:::.simulate,args=cbind(frame[i,],p=500))
#      foldid <- rep(c(0,1),times=c(frame$n[i],10000))
#      set.seed(i)
#      ridge[[i]] <- do.call(what=cornet:::cv.cornet,args=c(data,alpha=0,foldid=foldid))
#      set.seed(i)
#      lasso[[i]] <- do.call(what=cornet:::cv.cornet,args=c(data,alpha=1,foldid=foldid))
#  }
#  names(lasso) <- names(ridge) <- paste0("set",seq_len(nrow(frame)))
#  save(lasso,ridge,frame,file="results/simulation.RData")
#  
#  writeLines(text=capture.output(utils::sessionInfo(),cat("\n"),
#          sessioninfo::session_info()),con="results/info_sim.txt")

## ----figure_BOX,eval=FALSE----------------------------------------------------
#  #--- boxplot of different metrics ---
#  load("results/simulation.RData",verbose=TRUE)
#  
#  fuse0 <- fuse1 <- list()
#  for(i in c("deviance","class","mse","mae","auc")){
#    fuse0[[i]] <- sapply(ridge,function(x) (x[[i]]["combined"]-x[[i]]["binomial"]))
#    fuse1[[i]] <- sapply(lasso,function(x) (x[[i]]["combined"]-x[[i]]["binomial"]))
#  }
#  
#  grDevices::pdf("manuscript/figure_BOX.pdf",width=6,height=4)
#  graphics::par(mar=c(1.9,1.9,0.1,0.1))
#  graphics::plot.new()
#  ylim <- range(unlist(fuse0),unlist(fuse1))
#  at <- seq(from=1,to=9,by=2)
#  graphics::plot.window(xlim=c(min(at)-0.6,max(at)+0.6),ylim=ylim)
#  graphics::axis(side=2)
#  graphics::abline(h=0,col="grey",lty=2)
#  graphics::abline(v=at+1,col="grey",lty=2)
#  graphics::box()
#  graphics::boxplot(fuse1,at=at-0.5,add=TRUE,axes=FALSE,col="white",border="black")
#  graphics::boxplot(fuse0,at=at+0.5,add=TRUE,axes=FALSE,col="white",border="darkgrey")
#  labels <- names(fuse1)
#  labels <- ifelse(labels=="class","mcr",labels)
#  labels <- ifelse(labels %in% c("mcr","mse","mae","auc"),toupper(labels),labels)
#  for(i in seq_along(labels)){
#    graphics::axis(side=1,at=at[i],labels=bquote(Delta ~ .(labels[i])))
#  }
#  grDevices::dev.off()
#  
#  # decrease
#  sapply(fuse1,function(x) mean(x<0)) # lasso
#  sapply(fuse0,function(x) mean(x<0)) # ridge
#  
#  # constant
#  sapply(fuse1,function(x) mean(x==0)) # lasso
#  sapply(fuse0,function(x) mean(x==0)) # ridge
#  
#  # increase
#  sapply(fuse1,function(x) mean(x>0)) # lasso
#  sapply(fuse0,function(x) mean(x>0)) # ridge

## ----figure_TAB,eval=FALSE----------------------------------------------------
#  #--- plot of percentage changes ---
#  load("results/simulation.RData",verbose=TRUE)
#  
#  loss <- list()
#  loss$ridge <- as.data.frame(t(sapply(ridge,function(x) x$deviance)))
#  loss$lasso <- as.data.frame(t(sapply(lasso,function(x) x$deviance)))
#  
#  data <- list()
#  for(i in c("ridge","lasso")){
#    data[[i]] <- data.frame(row.names=rownames(frame))
#    data[[i]]$"(1)" <- 100*(loss[[i]]$binomial-loss[[i]]$intercept)/loss[[i]]$intercept
#    data[[i]]$"(2)" <- 100*(loss[[i]]$combined-loss[[i]]$intercept)/loss[[i]]$intercept
#    data[[i]]$"(3)" <- 100*(loss[[i]]$combined-loss[[i]]$binomial)/loss[[i]]$binomial
#  }
#  
#  row <- colnames(data$lasso)
#  col <- colnames(frame)
#  txt <- expression(rho,n,s,sigma,t,q)
#  
#  for(k in c("ridge","lasso")){
#    grDevices::pdf(paste0("manuscript/figure_",k,".pdf"),width=6.5,height=4)
#    graphics::par(mfrow=c(length(row),length(col)),
#                mar=c(0.2,0.2,0.2,0.2),oma=c(4,4,0,0))
#    for(i in seq_along(row)){
#      for(j in seq_along(col)){
#        y <- data[[k]][[row[i]]]
#        x <- frame[[col[j]]]
#        graphics::plot.new()
#        graphics::plot.window(xlim=range(x),ylim=range(y),xaxs="i")
#        graphics::box()
#        graphics::abline(h=0,lty=1,col="grey")
#        graphics::points(y=y,x=x,cex=0.5,pch=16,col=ifelse(y>0,"black","grey"))
#        line <- stats::loess.smooth(y=y,x=x,evaluation=200)
#        graphics::lines(x=line$x,y=line$y,col="black",lty=2,lwd=1)
#        if(j==1){
#          graphics::mtext(text=row[i],side=2,line=2.5,las=2)
#          graphics::axis(side=2)
#        }
#        if(i==length(row)){
#          graphics::mtext(text=txt[j],side=1,line=2.5)
#          graphics::axis(side=1)
#        }
#      }
#    }
#    grDevices::dev.off()
#  }
#  
#  cbind(col,as.character(txt)) # verify

