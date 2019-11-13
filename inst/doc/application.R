## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo=TRUE,eval=FALSE)
#setwd("C:/Users/armin.rauschenberger/Desktop/cornet") # local drive
#setwd("Z:/Rauschenberger/cornet/BIOINF_20XX-XX-XX") # shared drive
#utils::install.packages(pkgs=c("devtools","missRanger","xtable"))
#devtools::install_github("rauschenberger/cornet")

## ----process-------------------------------------------------------------
#  # features
#  X <- read.csv("data/PPMI_Baseline_Data_02Jul2018.csv",row.names="PATNO",na.strings=c(".",""))
#  X <- X[X$APPRDX==1,] # Parkinson's disease
#  X[c("SITE","APPRDX","EVENT_ID","symptom5_comment")] <- NULL
#  100*mean(is.na(X)) # proportion missingness
#  #x <- mice::complete(data=mice::mice(X,m=10,maxit=5,method="pmm",seed=1),action="all") # low-dimensional
#  x <- lapply(seq_len(10),function(x) missRanger::missRanger(data=X,pmm.k=3,
#          num.trees=100,verbose=0,seed=1)) # high-dimensional
#  x <- lapply(x,function(x) model.matrix(~.-1,data=x))
#  
#  # outcome
#  Y <- read.csv("data/PPMI_Year_1-3_Data_02Jul2018.csv",na.strings=".")
#  Y <- Y[Y$APPRDX==1 & Y$EVENT_ID %in% c("V04","V06","V08"),]
#  Y <- Y[,c("EVENT_ID","PATNO","moca")]
#  Y <- reshape(Y,idvar="PATNO",timevar="EVENT_ID",direction="wide")
#  rownames(Y) <- Y$PATNO; Y$PATNO <- NULL
#  
#  # overlap
#  names <- intersect(rownames(X),rownames(Y))
#  Y <- Y[names,]; x <- lapply(x,function(x) x[names,])
#  save(Y,x,file="data/processed_data.RData")
#  writeLines(text=capture.output(utils::sessionInfo(),cat("\n"),
#          sessioninfo::session_info()),con="data/info_data.txt")
#  
#  # dictionary
#  rm(list=ls())
#  info <- read.csv("data/PPMI_Data_Dictionary_for_Baseline_Dataset_Jul2018.csv",nrows=238)
#  info <- info[!info$Variable %in% c("","PATNO","SITE","APPRDX","EVENT_ID","symptom5_comment"),]
#  info <- paste(paste0("\\textit{",info$Variable,"}: ",info$Description),collapse="; ")
#  info <- gsub(pattern=" \\(OFF\\)",replacement="",x=info)
#  info <- gsub(pattern="_",replacement="\\_",x=info,fixed=TRUE)
#  info <- gsub(pattern="&",replacement="\\\\&",x=info)
#  cat(info)

## ----analyse-------------------------------------------------------------
#  load("data/processed_data.RData",verbose=TRUE)
#  
#  colSums(!is.na(Y)) # sample size
#  round(100*colMeans(Y<25.5,na.rm=TRUE),1) # proportion impairment
#  
#  # --- Cross-validating models. ---
#  names <- paste0(rep(c("lasso","ridge"),each=3),seq_len(3)) # change to 3!
#  loss <- fit <- pval <- list()
#  for(i in seq_along(x)){
#    loss[[i]] <- fit[[i]] <- pval[[i]] <- list()
#      cat("i =",i,"\n")
#      for(j in seq_along(names)){
#        cat("j =",j," ")
#        alpha <- 1*(substr(names[j],start=1,stop=5)=="lasso")
#        index <- as.numeric(substr(names[j],start=6,stop=6))
#        y <- Y[,index]
#        cond <- !is.na(y)
#        set.seed(i)
#        loss[[i]][[names[j]]] <- cornet::cv.cornet(y=y[cond],cutoff=25.5,
#                                    X=x[[i]][cond,],alpha=alpha)
#        set.seed(i)
#        fit[[i]][[names[j]]] <- cornet::cornet(y=y[cond],cutoff=25.5,
#                                    X=x[[i]][cond,],alpha=alpha)
#        set.seed(i)
#        pval[[i]][[names[j]]] <- median(replicate(n=50,expr=cornet:::.test(y=y[cond],cutoff=25.5,X=x[[i]][cond,],alpha=alpha)))
#      }
#      cat("\n")
#  }
#  
#  save(loss,fit,pval,file="results/application.RData")
#  writeLines(text=capture.output(utils::sessionInfo(),cat("\n"),
#          sessioninfo::session_info()),con="results/info_app.txt")

## ----compare-------------------------------------------------------------
#  load("results/application.RData",verbose=TRUE)
#  
#  names <- c(paste0("lasso",1:3),paste0("ridge",1:3))
#  frame <- data.frame(row.names=names)
#  for(i in names){
#  
#    # deviance
#    dev <- as.data.frame(t(sapply(loss,function(x) x[[i]]$deviance)))
#    dec <- (dev$combined-dev$binomial)/dev$binomial
#    # change in percent
#    frame[i,"dev.min"] <- round(100*min(dec),digits=1)
#    frame[i,"dev.max"] <- round(100*max(dec),digits=1)
#  
#    # proportion with improvement
#    frame[i,"dev.num"] <- sum(dev$combined < dev$binomial)
#  
#    # significance based on multi-split
#    pvalue <- sapply(pval,function(x) x[[i]])
#    frame[i,"pval.min"] <- signif(min(pvalue),digits=2)
#    frame[i,"pval.max"] <- signif(max(pvalue),digits=2)
#    #frame[i,"pval.prop"] <- round(mean(pvalue<0.05),digits=2)
#  
#    # quantiles weight parameter
#    q <- round(quantile(sapply(fit,function(x) x[[i]]$pi.min),probs=c(0,0.5,1)),digits=2)
#    frame[i,"pi.min"] <- q[1]
#    frame[i,"pi.med"] <- q[2]
#    frame[i,"pi.max"] <- q[3]
#  
#    # quantiles scale parameter
#    q <- round(quantile(sapply(fit,function(x) x[[i]]$sigma.min),probs=c(0,0.5,1)),digits=2)
#    frame[i,"sd.min"] <- q[1]
#    frame[i,"sd.med"] <- q[2]
#    frame[i,"sd.max"] <- q[3]
#  }
#  
#  # presentation
#  
#  rownames(frame) <- paste0(substr(x=rownames(frame),start=1,stop=5)," ",
#                            substr(x=rownames(frame),start=6,stop=6))
#  colnames(frame) <- gsub(pattern="dev.",replacement="\\\\delta_{\\\\text{",x=colnames(frame))
#  colnames(frame) <- gsub(pattern="pval.",replacement="p_{\\\\text{",x=colnames(frame))
#  colnames(frame) <- gsub(pattern="pi.",replacement="\\\\pi_{\\\\text{",x=colnames(frame))
#  colnames(frame) <- gsub(pattern="sd.",replacement="\\\\sigma_{\\\\text{",x=colnames(frame))
#  colnames(frame) <- paste0("$",colnames(frame),"}}$")
#  
#  digits <- 2+c(0,1*grepl(pattern="delta_|p_",x=colnames(frame))) # was "p_"
#  xtable <- xtable::xtable(frame,digits=digits)
#  xtable::print.xtable(xtable,include.rownames=TRUE,sanitize.text.function=function(x) x)

## ----figure_MAP----------------------------------------------------------
#  load("results/application.RData",verbose=TRUE)
#  sum <- fit[[1]]$lasso1
#  sum$cvm <- Reduce("+",lapply(fit,function(x) x$lasso1$cvm))
#  sum$sigma.min <- sapply(fit,function(x) x$lasso1$sigma.min)
#  sum$pi.min <- sapply(fit,function(x) x$lasso1$pi.min)
#  
#  grDevices::pdf("manuscript/figure_MAP.pdf",width=4,height=3)
#  #grDevices::postscript("manuscript/figure_MAP.eps",horizontal=FALSE,onefile=FALSE,paper="special",width=4,height=3)
#  graphics::par(mar=c(4,4,0.5,0.5))
#  cornet:::plot.cornet(sum)
#  grDevices::dev.off()

## ----figure_TFN----------------------------------------------------------
#  rm(list=ls())
#  #load("results/application.RData",verbose=TRUE)
#  #sigma.min <- sapply(fit,function(x) x$lasso1$sigma.min)
#  #sigma.min <- sapply(fit,function(x) sapply(x,function(x) x$sigma.min))
#  #sigma <- round(quantile(sigma.min,p=c(0,0.5,1)),digits=2)
#  #cutoff <- unique(sapply(fit,function(x) x$lasso1$cutoff))
#  
#  sigma <- c(1,2,3); cutoff <- 25.5 # or actual values (see above)
#  x <- seq(from=20,to=30,length.out=100)
#  
#  grDevices::pdf("manuscript/figure_TFN.pdf",width=4,height=3)
#  #grDevices::postscript("manuscript/figure_TFN.eps",horizontal=FALSE,onefile=FALSE,paper="special",width=4,height=3)
#  graphics::par(mar=c(4,4,0.5,0.5))
#  graphics::plot.new()
#  graphics::plot.window(xlim=range(x),ylim=c(0,1))
#  graphics::box()
#  graphics::axis(side=1)
#  graphics::axis(side=2)
#  graphics::title(xlab=expression(hat(y)),ylab=expression(Phi(hat(y),mu,sigma^2)),line=2.5)
#  graphics::abline(h=0.5,lty=2,col="grey")
#  graphics::abline(v=cutoff,lty=2,col="grey")
#  
#  lty <- c(2,1,3); lwd <- c(1,1,2)
#  lty <- c("dashed","solid","dotted")
#  for(i in seq_along(sigma)){
#    p <- stats::pnorm(q=x,mean=cutoff,sd=sigma[i])
#    graphics::lines(x=x,y=p,lty=lty[i],lwd=lwd[i])
#  }
#  
#  graphics::text(x=cutoff,y=0.40,labels=bquote(mu==.(cutoff)),pos=4)
#  legend <- sapply(sigma,function(x) as.expression(bquote(sigma == .(x))))
#  graphics::legend(x="topleft",legend=legend,lty=lty,bty="n",lwd=lwd)
#  grDevices::dev.off()
