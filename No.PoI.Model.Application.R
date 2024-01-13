#install.packages("readxl")
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
p=248
alpha=.05
grd <- seq(0,1,len=248)



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

grd2ind <- function(grd_value, grd){
  ## Projects a grid value of a point to its index
  apply(as.matrix(grd_value), 1, function(x) which.min(abs(grd-x)))
}
FunRegPoISim <- function(N, grd, DGP = DGP_name,
                         error_sd = sigma) {
  
  mean_Y   <- 0		# Mean of Y
  sd_Y     <- sigma        # Sd of scalar Y errors
  
  ## set the PoI-parameters (tau, PoI-coefs , beta function)
  PoIDGP     <- setPoIDGP(DGP_name)
  beta_fun   <- function(t){ eval(parse(text=PoIDGP$fct_text))}
  beta_fun_grd <- beta_fun(t=grd)
  tau_ind_true <- grd2ind(PoIDGP$tau, grd)    
  ## Evaluate the X functions at the true taus
  X_tau <- t(X_mat[tau_ind_true, 1:N, drop=F]) # N x length(tau)
  
  ## PoI Model: Y_i = beta_0 + int X_i(t) beta(t) + sum_j beta_j X_i(tau_j) + eps_i
  Y     <-  t(X_mat) %*% beta_fun_grd/p + X_tau %*% PoIDGP$beta_tau + rnorm(N, mean = mean_Y, sd = sd_Y)
  
  simObj <- list(Y = Y ,
                 X = X_mat,
                 trueBetaPoI = PoIDGP$beta_tau,
                 trueBetaCurve = beta_fun_grd,
                 trueTauGrd = PoIDGP$tau,
                 trueTauInd = tau_ind_true
  )
  class(simObj) <- "FunRegPoISim"
  simObj
}
ind2grd <-
  function(ind, grd){
    ## Projects a index of a point to a grid value
    return(grd[ind])
  }
calcScoreBetaFun <-
  function(ef, coefs , p){
    ## returns a function object from (evaluated) eigenfunctions and coefficients
    if ( length(coefs) == 1 ){
      y  <- ef * coefs
    } else{
      y  <- ef %*% coefs
    }
    return(splinefun( y =  y, x = seq(0,1 , len=p)))
  }
calcPoIandScoreLM <-
  function(estDataFrame, s, k, K, S){
    ## ########################
    ## A Function calculating for arbitrary s PoIs the corresponding lm models and their BICs
    ## Input: # k =1; s = 1;K = maxK; S = maxS
    ##  - s: the number of included PoIs
    ##	- k: the number of included eigenfunctions
    ##  - estDataFrame: a named data.frame with 1 + K + S columns of Y, the K scores and S potential PoIs
    ##  - K: the maximum number of Eigenfunctions
    ##  - S: the maximum number of PoIs
    ##  - model_out: should the lm fit be returned?
    ## Output:
    ##  - The MLR BIC
    dep_name  <- names( estDataFrame )[1]
    
    ## print(paste0("s = ", s, " k = ", k , " K = ", K, " S = ", S))
    if(s > 0 && s <= S && k <= K){
      ## Get the dependent variable s name
      ef_names   <- NULL 
      poi_names  <- NULL 
      if ( k > 0 ){ ## If other regressors are in estDataFrame, get the name of them
        ef_names <- (names( estDataFrame )[-1])[1:k]	
      } 
      if ( s > 0 ){ ## Get PoI names
        poi_names <- (names( estDataFrame )[-(1:(K+1))])[1:s]
      } 
      ## Generate a formula depending on the number of chosen PoI and estimate it via LM
      s_formula    <- as.formula(paste0( dep_name, " ~ -1 +" , paste( c(ef_names, poi_names) , collapse=" + ")))
      lmfit        <- lm(formula = s_formula , data = estDataFrame)
      BIC          <- calcBIClm(lmfit)[2]
      #GCV          <- calcGCVlm(lmfit)[2]
      return(list(BIC = BIC, #GCV = GCV,
                  lmFit = lmfit))
      
    } else if ( s == 0 ) {
      ef_names    <- NULL 
      if ( k > 0 ){ ## If other regressors are in estDataFrame, get the name of them
        ef_names  <- (names( estDataFrame )[-1])[1:k]
        s_formula <- as.formula(paste0( dep_name, " ~ -1+" , paste(ef_names , collapse = " + ")))
      } 
      else {
        s_formula <- as.formula(paste0( dep_name, " ~ -1" ))
      }
      lmfit  <- lm(formula = s_formula , data = estDataFrame)
      
      BIC    <- calcBIClm(lmfit)[2]
      #GCV    <- calcGCVlm(lmfit)[2]
      return(list(BIC = BIC, #GCV = GCV,
                  lmFit = lmfit))
      
    } else{
      stop("Something is wrong with estDataFrame with s and k")
    }
  }
calcBIClm <-
  function(fit){
    ## ###########################################
    ## A Function caluclating the BIC of an lm fit
    n   <- length(fit$residuals)
    edf <- n - fit$df.residual
    RSS <- sum(weighted.residuals(fit)^2)
    bic <- n * log(RSS/n) + edf * log(n)
    if (bic==-Inf) {
      stop("BIC is -Inf")
    } else {
      return(c(edf, bic))
    }
  }
createEstDFforDirectSearch <-
  function(Y, X, PoI_vals){
    ## ########################
    ## A Function creating the DF for the DirectSearch Function
    ## Input:
    ## -  Y
    ## -  X
    ## -  PoI_vals
    Y_names <- "Y"
    
    ## if no PoIs
    if( length(PoI_vals) == 0 ){ 	
      if( is.null(X) ){ 			# and if no regressors
        estDataFrame  <- data.frame(Y)
        return(estDataFrame)
      } else { 					# and if regressors
        X_names             <- paste("PC", 1:dim(X)[2], sep="")
        estDataFrame        <- data.frame( cbind(Y , X) )
        names(estDataFrame) <- c( Y_names , X_names  )
        return(estDataFrame)
      }
      ## if exactly one PoI
    } else if (length(PoI_vals) == length(Y) ){ 
      if( is.null(X) ){ 			# and if no regressors
        PoI_names           <- paste("PoI.1")
        estDataFrame        <- data.frame(cbind(Y , PoI_vals))
        names(estDataFrame) <- c( Y_names , PoI_names  )
        return(estDataFrame)
      } else {  					# and if regressors
        X_names             <- paste("PC", 1:dim(X)[2], sep="")
        PoI_names           <- paste("PoI.1")
        estDataFrame        <- data.frame(cbind(Y , X, PoI_vals))
        names(estDataFrame) <- c( Y_names , X_names, PoI_names  )
        return(estDataFrame)
      }
      ## if multiple PoIs 
    } else if (dim(PoI_vals)[1] == length(Y) ){ 
      if( is.null(X) ){ 			# and if no regressors
        PoI_names           <- paste("PoI", 1:dim(PoI_vals)[2], sep=".")
        estDataFrame        <- data.frame(cbind(Y , PoI_vals))
        names(estDataFrame) <- c( Y_names , PoI_names  )
        return(estDataFrame)
      } else {  					# and if regressors
        if(length(X) == length(Y)){ # if 1 regressors
          X_names <- "PC1"
        } else{
          X_names  <- paste("PC", 1:dim(X)[2], sep="")
        }
        PoI_names    <- paste("PoI", 1:dim(PoI_vals)[2], sep=".")
        estDataFrame <- data.frame(cbind(Y , X, PoI_vals))
        names(estDataFrame) <- c( Y_names , X_names, PoI_names  )
        return(estDataFrame)
      }	
    } else {
      stop("Something is wrong with the dimensions in X and Y")
    }
  }
searchPotPoi <-
  function(k, X, y, dom = c(0,1), plotting=FALSE){
    ## ########################
    ## A Function identifying potential points of impact of Cov(Z_delta, Y)
    ## Input:
    ##  - k: grid indices on domain
    ##  - X: centered X-data; Nxp MMatrix
    ##  - y: centered y-data
    ##  - dom: the domain where X(t) s live on
    ##  - plotting: binary variable if covariance should be plotted
    ## Output:
    ##  - ind.x: index of potential points of impact
    ##X <- X_st; y = Y ; k = 1
    N           <- ncol(X)
    p 		<- nrow(X)              	# number of discretization points 
    a 		<- range(dom)[1]			# left border
    b 		<- range(dom)[2]			# right border
    t 		<- seq(a,b,length.out=p) 	# grid in the domain [a,b]
    
    delta  	<- k*(b-a)/(p-1)    		# translate k to delta
    Jdelta 	<- (k+1):(p-k)      		# only 'inner' grid-point in [a,b]
    ## print(X)
    ## print(y)
    CorXY      <- 1/N*((X%*%y)) #up to factor 1/n
    CorZY      <- apply(CorXY,2,function(X){X[Jdelta]-1/2*( X[Jdelta-k] + X[Jdelta+k])})
    
    COR        <- abs(CorZY) #Criteria
    COR.aux    <- COR
    Jdelta.aux <- Jdelta
    t.delta    <- t[Jdelta]
    t.aux      <- t.delta
    taucand    <- vector() #container for hat(tau)
    ## Estimate tau!
    ##print(COR)
    while(sum(COR.aux>0.000000001)>1){
      tau.aux <- which.max(COR.aux)  
      taucand <- c(taucand,t.aux[tau.aux])
      indi    <- (abs(t.aux-t.aux[tau.aux])>=sqrt(delta)/2)# theoretical threshold
      t.aux   <- t.aux[indi]
      COR.aux <- COR[t.delta%in%t.aux]
    }
    if(plotting){
      plot(t.delta,COR,xlab="",ylab="",yaxt='n',type="l")#main=paste("k =",k, "kappa.est =",round(kappa,2), "   c =", c, "delta= ", round(delta,3)))
      ##abline(v=t[tau],col="black",lwd=1)
      abline(v=taucand,col="red")
    }
    ## Get the indexsets corresponding to the tau-candidates and return them:
    ind.x <- which(t%in%taucand)[rank(taucand)]
    
    list(ind.x = ind.x, cor = COR)
  }
KnePosSarEstimation_dirSearch <-
  function(Y, X_mat, add.vars, k, grd, dom, threshold_cpv = 0.95, plotting = FALSE, maxK = 40, maxPoI=0){
    ## ########################
    ## A Function performing the KnePosSar algorithm of doing the PoI esimation by directed BIC search
    ## Input:
    ##  - Y: scalar dependent variables (N)
    ##  - X_mat: data matrix (p x N)
    ##  - k: the number of grid indices to to use for potential PoISearch; k <- 1
    ##  - grd: the grid of length p where the functions are observed
    ##  - dom: the domain, where the function is observed
    ##  - plotting=FALSE: should be plots provided?
    ##  - threshold_cpv = FALSE
    ## Output:
    ##  - k: the k which was used
    ##  - selS: the chosen number of PoIs
    ##  - estPoI: the estimated Parameters of chosen PoIs
    ##  - estTau: the estimated taus on the domain
    ##  - estTauGrd: the estimated taus on the grid
    ##  - Y_hat: the predicted values
    ##  - estBeta: beta_hat(t) function
    ##  - BIC: the minimal BIC after which the model was chosen
    
    ## todo (additional variables currently only implemented in CraKneSaPoI):
    add.var.coef <- NULL
    
    ## Get potential points of impact
    p      <- nrow(X_mat)
    Y_st   <- scale(Y)
    X_st   <- scale(X_mat) 
    potPoI <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting=plotting) 
    
    ## Take minimum of up to maxPoI or potPoi number PoI into your choice
    maxS   <- min(length(potPoI$ind.x), maxPoI)
    if(maxS > 0){
      potPoI <- potPoI$ind.x[1:maxS]
    } else{
      potPoI <- NULL
    }
    
    ## Centering of Y and X data
    X_c	<- apply( X_mat, 1 , FUN = function(row) return(row - mean(row)) )
    Y_c	<- Y - mean(Y)
    
    ## Eigenfunctions, Basis Coefficients and Eigenvalues via PCA
    efSelection  <- selectEFwithCumPropVar(X_mat = X_mat, threshold_cpv = threshold_cpv)
    scores       <- X_c %*% efSelection$evecs * 1/p  
    
    ## If CPV is supplied, estimate the number of k via CPV using the selection by the selectEFwithCumPropVar function
    if( is.numeric(threshold_cpv) ) { # 
      ## If CPV is supplied, estimate the number of k via CPV
      K_choice <- efSelection$K
      
      ## Prepare Data for BIC Choice
      estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ ,1:K_choice], PoI_vals = X_c[ ,potPoI, drop=FALSE])
      
      # A directed S search over all number of S; choose model via minimal BIC
      LMandBICs <- lapply(0:maxS, function(s) calcPoIandScoreLM(estDataFrame = estDataFrame, s = s, k = K_choice, K = K_choice, S = maxS))
      S_choice  <- (0:maxS)[ which.min(sapply(LMandBICs, function(x) x[["BIC"]])) ] # S_choice = 5
      BICChoice <- min(sapply(LMandBICs, function(x) x[["BIC"]]))
      
      ## If no cpv is supplied, estimate the number of K by dirSearch	and the number S_choice by fullSearch
    } else if( isTRUE(!threshold_cpv) ) { 
      
      ## Prepare Data for directed BIC Search over all K and all S head(estDataFrame)
      estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ ,1:maxK], PoI_vals = X_c[ ,potPoI, drop=FALSE])
      
      ## A directed S search over all number of S and all number of K up to maxK
      ## For each K generate a list of all BICs per s in S
      
      ## BICs <- lapply(0:maxK, function(k) { sapply(0:maxS, function(s) calcPoIandScoreLM(estDataFrame = estDataFrame, s = s, k = k, K = maxK, S = maxS)$BIC)})
      ## way faster:
      BICs <- lapply(0:maxK , function(k){
        vapply( 0:maxS, function(s) {
          if (k==0){
            if (s==0){
              y <- estDataFrame[,1]
              n <- length(y)
              RSS <- sum(y^2)
              return(n * log(RSS/n))  # here: edf = 0
            } else {
              lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c((maxK+2):(maxK+s+1))]))
            }
          } else {
            if (s==0) {
              lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, 2:(k+1) ]))
            } else {
              lmEst <- lm.fit(y = estDataFrame[,1], x = as.matrix(estDataFrame[, c(2:(k+1),(maxK+2):(maxK+s+1))]))
            }
          }
          calcBIClm(lmEst)[2]
        }, 0)
      })
      
      ## Generate for each k optimal included 
      BICs_per_k <- cbind(K_ind = 0:maxK, t(sapply(BICs, function(bic_k) return(list(S_ind = which.min(bic_k), K_BIC = min(bic_k))))))
      optSoptK   <- BICs_per_k[  which.min(BICs_per_k[ , "K_BIC"]) , c("K_ind", "S_ind", "K_BIC")]
      BICChoice  <- unlist(optSoptK)["K_BIC"]
      S_choice   <- (0:maxS)[ optSoptK[["S_ind"]] ] 		# If optSoptK["S_ind"]=1 --> (0:maxS)[ (optSoptK[["S_ind"]]) ] = 0 -> adjustment
      K_choice   <- (0:maxK)[ (optSoptK[["K_ind"]]+1) ] 	# K_choice = 3
    } else {
      stop("Problem using CPV")
    }
    
    if(S_choice > 0){
      PoIChoice <- potPoI[1:S_choice]
    } else {
      PoIChoice <- NULL
    }
    
    ## get the final model estimates names(estDataFrame)
    if(K_choice != 0){
      estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = scores[ , 1:K_choice, drop=FALSE], PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    } else {
      estDataFrame <- createEstDFforDirectSearch(Y = Y_c, X = NULL, PoI_vals = X_c[ , PoIChoice, drop=FALSE])
    }
    
    ## estimate Y againgst eigenfunctions and PoI choices
    lm.fit     <- calcPoIandScoreLM(estDataFrame, s = S_choice, k = K_choice, K = K_choice, S = maxS)$lmFit
    coefnames  <- names(lm.fit$coeffic)
    score_names  <- coefnames[substr(coefnames , 1,3)!="PoI"]
    poi_names  <- coefnames[substr(coefnames , 1,3)=="PoI"]
    estEF      <- efSelection$evecs[ , score_names]
    estPsi     <- lm.fit$coefficients[score_names]
    estPoI     <- lm.fit$coefficients[poi_names]
    
    ## Calculate beta_hat(t), PoI_hat, Y_hat
    beta_hat   <- calcScoreBetaFun(ef = estEF, coefs = estPsi, p=p)
    if (S_choice > 0) {
      ## Y_i = integral X_i(t) beta_hat(t) dt + sum_k beta_hat_k * X_i(tau_hat_k)
      Y_hat  <- (X_c %*% beta_hat(grd))/ p + X_c[, PoIChoice ,drop=FALSE] %*% estPoI
      
    } else {
      Y_hat  <- (X_c %*% beta_hat(grd))/ p
    }
    hat_matrix <- NULL
    list(k  = k, 
         selS = S_choice, 
         K_choice = K_choice,
         estPoI = estPoI,
         betaAddVar = add.var.coef,
         estTau = ind2grd(PoIChoice , grd),
         estTauGrd = PoIChoice,
         estPsi = estPsi,
         scores = scores ,
         estEF = estEF,
         Y_hat = Y_hat, 
         estBeta = beta_hat(grd),
         BIC = BICChoice,
         hat_matrix = hat_matrix
    )
  }

KnePosSarEstimation <-
  function(Y, X_mat, grd, add.vars, maxPoI = 0 ,threshold_cpv=0.95, estOrder = NULL, searchMethod = "dirSearch" ,
           k_seq = 1:floor(length(grd)/6) , maxK = 40, exPost = FALSE, scaleSearchPoI = TRUE,...){
    ## Function calculating the eigenfunction-based estimate of FunRegPoI. If estOrder ==NULL, the standard
    ## KPS-Estimation is done, in this case, exPost has no effect
    ##
    ##
    ## tocheck: effect of dom in KnePosSarEstimation, adjust in KnePosSarEstimation_R{1,2}!
    ## centering??
    ##
    ## if estOrder =="R2" and exPost =="FALSE": "PES"-estimation (with EF instead of Splines)
    ## if estOrder =="R2" and exPost =="TRUE": "PESES"-estimation (with EF instead of Splines)
    if (is.null(estOrder)){
      ## wrapper around KnePosSarEstimation_fullSearch and KnePosSarEstimation_dirSearch, already over all k
      dom <- range(grd)
      ## for each k in k_seq calculate the KPS
      KnePosSarcomResults <- lapply(k_seq, function(k) {
        if (searchMethod=="dirSearch") {
          KnePosSarEstimation_dirSearch(Y = Y, X_mat = X_mat, add.vars = add.vars,k = k, grd = grd, dom = dom,
                                        threshold_cpv = threshold_cpv, plotting=FALSE, maxK = maxK )
        } else {
          KnePosSarEstimation_fullSearch(Y = Y, X_mat = X_mat, add.vars = add.vars, k = k, grd = grd, dom = dom,
                                         threshold_cpv = threshold_cpv, plotting=FALSE, maxK = maxK)
        }
      })
      
      bics <- lapply(KnePosSarcomResults , function(x) x$BIC)
      opt  <- which.min(unlist(bics))
      optEntry <- KnePosSarcomResults[[opt]]
      
    } else {
      
      Y_st   <- if (scaleSearchPoI) scale(Y) else Y - mean(Y)
      X_st  <- if (scaleSearchPoI) t(scale(t(X_mat))) else   t(scale(t(X_mat), scale = FALSE)) # dim(X_st)
      possible_deltas  <- ind2grd(k_seq, grd)
      
      ## get empirical eigenfunctions:
      efSelection  <- selectEFwithCumPropVar(X_mat = X_mat, threshold_cpv = threshold_cpv)
      
      
      ## for each k in k_seq calculate the KPS
      KnePosSarcomResults <- lapply(k_seq, function(k) {
        dom         <- range(grd)
        PoISearch   <- searchPotPoi(k, X_st, Y_st, dom=dom, plotting = FALSE) # with standardized values; might be better
        potPoI      <- PoISearch[["ind.x"]]
        cor         <- PoISearch[["cor"]]
        maxS        <- length(potPoI)
        
        if (estOrder=="R1"){
          KPSest <- KnePosSarEstimation_R1(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=potPoI,
                                           searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                           maxK = 40, maxPoI=0, efSelection = efSelection)
        } else {
          KPSest <- KnePosSarEstimation_R2(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=potPoI,
                                           searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                           maxK = 40, maxPoI=0, efSelection = efSelection)
        }
        KPSest[["cor"]]      <- PoISearch[["cor"]]
        KPSest[["k"]]        <- k
        KPSest[["delta"]]    <- possible_deltas[which(k == k_seq)]
        KPSest
      })
      names(KnePosSarcomResults)       <- possible_deltas
      
      if (exPost){
        exPostEstimation <- lapply(KnePosSarcomResults, function(entry) {
          
          if (estOrder=="R1"){
            KPSest <- KnePosSarEstimation_R1(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=entry[["estTauGrd"]],
                                             searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                             maxK = 40, maxPoI=0, efSelection = efSelection)
          } else {
            KPSest <- KnePosSarEstimation_R2(Y = Y, X_mat = X_mat, add.vars = add.vars, potPoI=entry[["estTauGrd"]],
                                             searchMethod = searchMethod, grd = grd, dom = dom, plotting = FALSE, 
                                             maxK = 40, maxPoI=0, efSelection = efSelection)
          }
          KPSest
        })
      }
      
      ## nest ExPost results in KnePosSarcomResults
      if (exPost){
        for( delta in possible_deltas){ # delta <- possible_deltas[1]
          delta_c <- which(delta == possible_deltas)
          KnePosSarcomResults[[delta_c]][["exPostEstimation"]] <- exPostEstimation[[delta_c]]
        }
      }
      ## Get Optima via BIC
      optimum   <- which.min(sapply(KnePosSarcomResults, function(x) x$BIC))
      
      if (exPost){
        exPostOptimum  <- which.min(sapply(KnePosSarcomResults, function(x) x[["exPostEstimation"]]$BIC))
      }
      
      ## save the selected entries
      KnePosSarcomResults[["opt"]]        <- list(KnePosSarcomResults[[optimum]], X = X_mat,Y = Y)
      names(KnePosSarcomResults[["opt"]]) <- c("res", "X", "Y")
      
      if (exPost){
        KnePosSarcomResults[["exPostOptimum"]] <-
          KnePosSarcomResults[[exPostOptimum]][["exPostEstimation"]]
      }
      ## get the selected estimation result:
      optEntry <- if (exPost) {
        KnePosSarcomResults[["exPostOptimum"]]
      } else {
        KnePosSarcomResults[["opt"]]$res
      }
      
    }
    
    
    
    if (!is.null(add.vars)){
      names(KnePosSarcomResults[[opt]]$betaAddVar) <-
        if (!is.null(names(add.vars))) {
          colnames(add.vars)
        } else {
          paste0("add.var", 1:col(add.vars))
        }
    }
    
    ## collect and make the return object
    estObj <- list(
      coefficients = list(
        betaCurve = as.vector( optEntry$estBeta) ,
        betaPoI = optEntry$estPoI ,
        betaAddVar = optEntry$betaAddVar,
        tauGrd = optEntry$estTau ,
        tauInd = optEntry$estTauGrd,
        selS = optEntry$selS) ,
      residuals = Y - as.vector(optEntry$Y_hat) - mean(Y) ,
      fitted = as.vector(optEntry$Y_hat) + mean(Y),
      data = list(Y = Y, X_mat = X_mat),
      call = "KPS" ,
      model = list(
        k = optEntry$k ,
        delta = ind2grd(optEntry$k , grd) ,
        cpv = threshold_cpv,
        maxPoI = maxPoI,
        K =  optEntry$K_choice,
        K_space= if(!missing(threshold_cpv)) NULL else c(1,maxK),
        estEigenFct = optEntry$estEF ,
        estScores = optEntry$scores,
        eff_df = optEntry$K_choice + optEntry$selS + length(add.vars) ,
        bic = optEntry$BIC , 
        bic_sq= NULL, ##
        hat_matrix= optEntry$hat_matrix,
        center = TRUE),
      kSeqRes = KnePosSarcomResults
    )
    return(estObj)
  }
FunRegPoI <-
  function(Y, X_mat, grd, add.vars = NULL, maxPoI= 0, ...) {
    ### Functional Linear Regression with Points of Impact
    ### wrapper around different estimators: PESES, KPS , CKS
    ###
    ### returns object of S3-class FunRegPoI, with methods plot and summary
    
    Y <- as.vector(Y)
    
    ## check for correct input
    if (!is.matrix(X_mat)) stop("'X_mat' has to be a matrix!")
    if (!is.null(add.vars) && !is.matrix(add.vars) && !is.vector(add.vars))
      stop("'add.vars' has to be a matrix or vector!")
    if (length(Y) != dim(X_mat)[2])
      stop("dimension of 'X_mat' and length of 'Y' don't match!")
    if (abs(diff(range(diff(grd)))) > 1e-10)
      stop("'grd' does not appear to be a regularly spaced grid!")
    
    ## if grid does not span 0-1, we scale X_mat accordingly
    if (abs(grd[1] - 1 / (2 * length(grd))) > 1e-10 ||
        abs(grd[length(grd)] + 1 / (2 * length(grd)) - 1) > 1e-10) {
      
      grid_scale <- (grd[length(grd)] - grd[1]) * (1 + 1 / length(grd))
      X_mat <- grid_scale * X_mat
    }
    
    estObj <- KnePosSarEstimation(Y, X_mat, grd, add.vars, maxPoI, ...)
    class(estObj) <- "FunRegPoI"
    estObj
  }

Est_KPS <- FunRegPoI(Y = Y , X_mat = X_mat , grd = grd , estimator = "KPS" , threshold_cpv = FALSE)
N <- ncol(X_mat)
Y_hat <- matrix(Est_KPS[3]$fitted,nrow=N)
prediction_error <- t(Y_hat-Y)%*%(Y_hat-Y)



