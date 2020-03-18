## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("cornet")

## ----eval=FALSE----------------------------------------------------------
#  #install.packages("devtools")
#  devtools::install_github("rauschenberger/cornet")

## ------------------------------------------------------------------------
library(cornet)

## ----eval=FALSE----------------------------------------------------------
#  set.seed(1)
#  n <- 100; p <- 500
#  X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#  beta <- rbinom(n=p,size=1,prob=0.05)
#  y <- rnorm(n=n,mean=X%*%beta)

## ----eval=FALSE----------------------------------------------------------
#  model <- cornet(y=y,cutoff=0,X=X)
#  model

## ----eval=FALSE----------------------------------------------------------
#  coef <- coef(model)

## ----eval=FALSE----------------------------------------------------------
#  predict <- predict(model,newx=X)

## ----eval=FALSE----------------------------------------------------------
#  cv.cornet(y=y,cutoff=0,X=X)

## ----eval=FALSE----------------------------------------------------------
#  #install.packages("BiocManager")
#  #BiocManager::install(c("GEOquery","Biobase"))
#  data <- GEOquery::getGEO(GEO="GSE80599")[[1]]
#  pheno <- Biobase::pData(data)
#  y <- as.numeric(pheno$`updrs-mds3.12 score:ch1`)
#  age <- as.numeric(pheno$`age at examination (years):ch1`)
#  gender <- ifelse(pheno$`gender:ch1`=="Female",1,0)
#  X <- cbind(age,gender,t(Biobase::exprs(data)))

## ----eval=FALSE----------------------------------------------------------
#  pvalue <- apply(X,2,function(x) cor.test(x,y)$p.value)
#  min(p.adjust(pvalue))
#  hist(pvalue)

## ----eval=FALSE----------------------------------------------------------
#  cor <- abs(cor(y,X,method="spearman"))
#  X <- X[,cor>0.3] # forbidden!

## ----eval=FALSE----------------------------------------------------------
#  #install.packages("BiocManager")
#  #BiocManager::install("GEOquery")
#  files <- GEOquery::getGEOSuppFiles("GSE97644")
#  pheno <- read.csv(textConnection(readLines(rownames(files)[1])))
#  y <- pheno$MOCA.Score
#  gender <- ifelse(pheno$Gender=="Female",1,0)
#  age <- pheno$Age
#  geno <- t(read.csv(textConnection(readLines(rownames(files)[2])),row.names=1))
#  X <- cbind(gender,age,geno)

## ----eval=FALSE----------------------------------------------------------
#  net <- cornet::cornet(y=y,cutoff=25,X=X)
#  set.seed(1)
#  cornet:::cv.cornet(y=y,cutoff=25,X=X)

## ----eval=FALSE----------------------------------------------------------
#  files <- GEOquery::getGEOSuppFiles("GSE95640")
#  X <- t(read.csv(textConnection(readLines(rownames(files)[1])),row.names=1))
#  y <- GEOquery::getGEO(GEO="GSE95640")[[1]] # no numeric outcome

## ----eval=FALSE----------------------------------------------------------
#  data <- GEOquery::getGEO(GEO="GSE109597")[[1]]
#  y <- as.numeric(Biobase::pData(data)$"bmi:ch1")
#  X <- t(Biobase::exprs(data))
#  cornet:::cv.cornet(y=y,cutoff=25,X=X,alpha=0)

## ----eval=FALSE----------------------------------------------------------
#  #install.packages("BiocManager")
#  #BiocManager::install("mixOmics")
#  set.seed(1)
#  data(liver.toxicity,package="mixOmics")
#  X <- as.matrix(liver.toxicity$gene)
#  Y <- liver.toxicity$clinic
#  cornet <- cornet::cornet(y=Y$BUN.mg.dL.,cutoff=15,X=X)
#  cornet:::cv.cornet(y=Y$BUN.mg.dL.,cutoff=15,X=X)
#  
#  loss <- list()
#  for(i in seq_along(Y)){
#    loss[[i]] <- cornet:::cv.cornet(y=Y[[i]],cutoff=median(Y[[i]]),alpha=0,X=X)
#  }
#  sapply(loss,function(x) x$deviance)

