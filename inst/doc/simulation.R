## ----setup,include=FALSE------------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,eval=FALSE)
#setwd("C:/Users/armin.rauschenberger/Desktop/cornet") # local drive
#setwd("Z:/Rauschenberger/cornet") # shared drive
#devtools::install_github("rauschenberger/cornet")

## ----analysis-----------------------------------------------------------------
#  iter <- 1000
#  set.seed(1)
#  frame <- data.frame(cor=runif(n=iter,min=0,max=0.9),
#                      n=round(runif(n=iter,min=100,max=200)),
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
#      set.seed(i)
#      ridge[[i]] <- do.call(what=cornet:::cv.cornet,args=c(data,alpha=0))
#      set.seed(i)
#      lasso[[i]] <- do.call(what=cornet:::cv.cornet,args=c(data,alpha=1))
#  }
#  names(lasso) <- names(ridge) <- paste0("set",seq_len(nrow(frame)))
#  save(lasso,ridge,frame,file="results/simulation.RData")
#  
#  writeLines(text=capture.output(utils::sessionInfo(),cat("\n"),
#          sessioninfo::session_info()),con="results/info_sim.txt")

## ----figure_BOX---------------------------------------------------------------
#  #--- boxplot of different metrics ---
#  load("results/simulation.RData",verbose=TRUE)
#  
#  fuse0 <- fuse1 <- list()
#  for(i in c("deviance","class","mse","mae","auc")){
#    fuse0[[i]] <- sapply(ridge,function(x) (x[[i]]["combined"]-x[[i]]["binomial"]))
#    fuse1[[i]] <- sapply(lasso,function(x) (x[[i]]["combined"]-x[[i]]["binomial"]))
#    #if(i=="auc"){fuse0[[i]] <- -fuse0[[i]]; fuse1[[i]] <- -fuse1[[i]]}
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
#  graphics::axis(side=1,at=at,labels=labels)
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

## ----figure_TAB---------------------------------------------------------------
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

