#install.packages("readxl")
#install.packages("leaps")
#install.packages("usethis")
#install.packages("devtools")
#install.packages("dplyr")
#install.packages("splines")
library("leaps")
library("usethis")
library("devtools")
library("dplyr")
library("splines")
library(readxl)
 X1 <- read_excel("D:/Temperature.xls")
 Y1 <- read_excel("D:/Humidity.xls")
X1 <- X1[,-c(1,2)]
X1 <- as.matrix(X1)
X1 <- as.numeric(X1)
X1 <- matrix(X1, nrow=248)
for(j in 1:ncol(X1)){
for(i in 1:nrow(X1)){
if(is.na(X1[i,j])==T){
X1[i,j] <- X1[i-1,j]
}
}
}
for(j in 1:ncol(X1)){
X1[,j] <- rev(X1[,j])
}
Y1 <- Y1[,-c(1,2)]
Y1 <- as.matrix(Y1)
for(j in 1:ncol(Y1)){
for(i in 1:nrow(Y1)){
if(is.na(Y1[i,j])==T){
Y1[i,j] <- Y1[i-1,j]
}
}
}
Y = apply(Y1,2,mean)
X_mat <- X1
X_mat <- X_mat-rowMeans(X_mat)
Y <- Y-mean(Y)
Y <- Y/100
devtools::load_all('C:\\Users\\Documents\\R\\win-library\\4.1\\FunRegPoI-master')
p=248
alpha=.05
grd <- seq(0,1,len=248)
Est_KPS <- FunRegPoI(Y = Y , X_mat = X_mat , grd = grd , estimator = "KPS" , threshold_cpv = FALSE)
(Est_KPS.sum <- summary(Est_KPS))
#plot(Est_KPS.sum)
plot(Est_KPS,xlab="Estimated points of impact",ylab="Estimated beta function")

selectEFwithCumPropVar <-
  function(X_mat, threshold_cpv = 0.95){
    p   <- nrow(X_mat)
    pca	<- prcomp(t(X_mat) )
    ## Adjust Eigenvectors/Eigenfunctions
    evecs <- pca$rotation * sqrt(p) # adjustment (approx) in order to have length 1 in L^2
    
    ## To check for the length of Eigenvectors
    ## apply(evecs, 1, FUN = function(row) return(sum((row-mean(row))^2)/p))
    
    ## Calc eigenvalues
    evals <- pca$sdev^2
    
    cpv   <- sapply(1:length(evals), FUN=function(k) return( sum(evals[1:k])/sum(evals)) )
    if( is.numeric(threshold_cpv) ) {
      K  <- which.max( cpv > threshold_cpv) 
    } else if( isTRUE(!threshold_cpv) ){
      K  <- dim(evecs)[2]
    }
    list( K = K, evals = evals, evecs = evecs, chosenEvecs = evecs[ ,1:K])
  }

k_sum       <- Est_KPS$model$k


index_tau    <-Est_KPS$coefficients$tauInd

Princi<-selectEFwithCumPropVar(X_mat = X_mat)

thethatij<-(t(X_mat)%*%Princi$evecs[,1:k_sum]) 
if((nrow(as.matrix(Princi$evecs[index_tau, 1:k_sum])))==(k_sum)){
  sumjkthetaphi<-thethatij%*%Princi$evecs[index_tau, 1:k_sum]/p
} else
{
  sumjkthetaphi<-thethatij%*%t(Princi$evecs[index_tau, 1:k_sum])/p
}

X_k_hat_tau<-X_mat[index_tau,]-t(sumjkthetaphi)
S_hat <- Est_KPS[1]$coefficients$selS
M<-matrix(rep(0, (S_hat)^2), nrow=(S_hat))
N <- ncol(X_mat)
for(i in 1:N){
  M<-M+X_k_hat_tau[,i]%*%t(X_k_hat_tau[,i])
}
beta_hat <- Est_KPS[1]$coefficients$betaPoI
numerator <- N*t(beta_hat)%*%M%*%beta_hat

Y_hat <- matrix(Est_KPS[3]$fitted,nrow=N)
denominator <- t(Y_hat-Y)%*%(Y_hat-Y)
T <- numerator/denominator 
H <- T-qgamma(1-(alpha), (S_hat)/2, 0.5, lower.tail = TRUE, log.p = FALSE)
print(H)

