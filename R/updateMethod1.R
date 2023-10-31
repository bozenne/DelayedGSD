## * Method1 (code)
updateMethod1 <- function(rho_alpha=2,          # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
                          rho_beta=2,           # rho parameter of the rho-family spending functions (Kim-DeMets) for beta
                          alpha=0.025, alphaSpent = NULL,         # Type-I error (overall) or spent up to the interim analysis (i.e. cumulative)
                          beta=0.2, betaSpent = NULL,            # Type-II error (overall) or spent at to the interim analysis (i.e. cumulative)
                          Kmax,                 # number of planned analyses (including the final analysis)
                          Info.max=NULL,        # Info.max, i.e. maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
                          InfoR.i=NULL,         # Expected or observed (wherever possible) information rates at the interim analyses 1:(Kmax-1)
                          InfoR.d=NULL,         # Expected or observed information rates at all potential decision analyses, including the final analysis 1:Kmax
                          uk = NULL,            # efficacy boundaries from the previous interim analyses or planning
                          lk = NULL,            # futility boundaries from the previous interim analyses or planning
                          k = NULL, type.k = NULL, ImaxAnticipated = FALSE, # current stage, type of analysis, and conclusion for all previous analyses
                          delta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
                          abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                          alternative="greater",   # greater is for Ho= theta > 0, "less" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider less)
                          binding=TRUE,         # whether the futility boundary is binding
                          Trace=FALSE,          # Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
                          nWhileMax=30,         # Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
                          toldiff= 1e-05,       # Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
                          tolcoef= 1e-04,       # Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
                          mycoefMax= 1.2,       # Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
                          mycoefL=1,            # Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
                          myseed=2902,           # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
                          cMin=-Inf,           # minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
                          PowerCorrection=FALSE #whether or not to apply a correction to the type II error spending to reduce the extent to which Method 1 is overpowered   
){
  ## {{{ set seed
  
    if(!is.null(myseed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }
        set.seed(myseed)
    }

  ## *** initialize boundaries
  
  #information sequence relevant for alpha spending and covariance matrix
  InfoR <- c(InfoR.i,InfoR.d[Kmax])
  
  ## *** compute variance-covariance matrix of vector (Z_1,...,Z_k)
  sigmaZk <- diag(1,Kmax)
  for(i in 1:Kmax){
    for(j in i:Kmax){
      sigmaZk[i,j] <- sqrt(InfoR[i]/InfoR[j])
      sigmaZk[j,i] <- sqrt(InfoR[i]/InfoR[j])
    }
  }
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  Info <- InfoR*Info.max
  
  #browser()
  
  # compute the mean of the multivariate normal distribution under the alternative H1
  thetheta <- delta*sqrt(Info)
  ## }}} 
  ## ** case k=1
    if(k==1){
        
    ## {{{ case k=1
    if(type.k=="interim"){
      alphaSpent[1] <- ErrorSpend(I=Info[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
      betaSpent[1] <-  ErrorSpend(I=Info[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
      
      uk[1] <- qnorm(p=1-alphaSpent[1],mean=0,sd=1)         # compute under the null (Ho)
      lk[1] <- qnorm(p=betaSpent[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)
    }
    
        ck.unrestricted <- calc_ck(uk=uk[1],
                                   lk=lk[1],
                                   Info.i=Info.i[1],
                                   Info.d=Info.d[1],
                                   Info.max=Info.max,
                                   ImaxAnticipated=ImaxAnticipated,
                                   rho_alpha=rho_alpha, 
                                   alpha=alpha,
                                   bindingFutility = binding)
                                        #------------
    }else{
    
    alphaSpentInc <- diff(c(0,alphaSpent))
    betaSpentInc <- diff(c(0,betaSpent))
    
    if(binding){
      TheLowerValues <- lk[1:(k-1)]
    }else{
      TheLowerValues <- rep(-Inf,k-1)
    }
    
    if(type.k=="interim"){
      
      ## ** Estimate uk
      alphaSpent[k] <- ErrorSpend(I=Info[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
      alphaSpentInc[k] <- alphaSpent[k] - alphaSpent[(k-1)]   
      betaSpent[k] <- ErrorSpend(I=Info[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
      
      if(PowerCorrection){
        #real type II error spent at analysis k-1
        betaSpentInc[k-1] <- TypeIIerrorSpent(lk=lk,
                                         uk=uk,
                                         ck=ck,
                                         Info.i=Info.i,
                                         Info.dk=Info.d[k-1],
                                         sigmaZk=sigmaZk,
                                         thetheta=thetheta,
                                         k=k-1,
                                         delta=delta,
                                         abseps=abseps)
        #correct the type II error that has been spent up to analysis k-1
        betaSpent[k-1] <- ifelse(k==2,betaSpentInc[k-1],betaSpent[k-2]+betaSpentInc[k-1]) 
      } 
      
      betaSpentInc[k] <- betaSpent[k] - betaSpent[(k-1)]   

        ## {{{ 
      ## {{{ u_k by solving what follows
      uk[k] <- uniroot(function(x){pmvnorm(lower = c(TheLowerValues,x),
                                           upper = c(uk[1:(k-1)],Inf),
                                           mean=rep(0,k),
                                           sigma= sigmaZk[1:k,1:k],
                                           abseps = abseps) - alphaSpentInc[k]},
                       lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)],  ## last boundary among the k-1 already computed that is not infinite  
                       upper = 1.1*uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)],
                       ## lower = lk[k-1],
                       ## upper = uk[k-1],
                       tol = abseps)$root
      
        ## ** Estimate lk
        lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
                                             upper = c(uk[1:(k-1)],x),
                                             mean=thetheta[1:k],
                                             sigma= sigmaZk[1:k,1:k],
                                             abseps = abseps) - betaSpentInc[k]},
                         lower = lk[k-1], 
                         upper = uk[k], 
                         tol = abseps)$root
    }
        if(type.k %in% c("interim","decision")){
            ck.unrestricted <- calc_ck(uk=uk[1:k],
                                       lk=lk[1:k],
                                       Info.i=Info.i[1:k],
                                       Info.d=Info.d[k],
                                       Info.max=Info.max,
                                       ImaxAnticipated=ImaxAnticipated,
                                       rho_alpha=rho_alpha,
                                       alpha=alpha,
                                       bindingFutility = binding)
        }else if(type.k=="final"){
            alphaSpent[k] <- alpha
            alphaSpentInc[k] <- alphaSpent[k] - alphaSpent[(k-1)]   
            betaSpent[k] <- beta ## ErrorSpend(I=Info.i[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
            betaSpentInc[k] <- betaSpent[k] - betaSpent[(k-1)]
            lowerRoot <- lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)]  ## last boundary among the k-1 already computed that is not infinite
            upperRoot <- 1.1*uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)]
            if(InfoR.d[k]>1){
                lowerRoot <- -10
                upperRoot <- 10
            }

            uk[k] <- uniroot(function(x){pmvnorm(lower = c(TheLowerValues,x),
                                           upper = c(uk[1:(k-1)],Inf),
                                           mean=rep(0,k),
                                           sigma= sigmaZk[1:k,1:k],
                                           abseps = abseps) - alphaSpentInc[k]},
                       lower = lowerRoot,
                       upper = upperRoot,
                       tol = abseps)$root
      
      ## lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
      ##                                      upper = c(uk[1:(k-1)],x),
      ##                                      mean=thetheta[1:k],
      ##                                      sigma= sigmaZk[1:k,1:k],
      ##                                      abseps = abseps) - betaSpentInc[k]},
      ##                  lower = lk[utils::tail(intersect(which(!is.infinite(lk)),1:(k-1)),1)],  ## last boundary among the k-1 already computed that is not infinite  
      ##                  upper = uk[utils::tail(intersect(which(!is.infinite(uk)),1:(k-1)),1)],  
      ##                  tol = abseps)$root
      
      lk[k] <- uk[k]
      
      ck.unrestricted <- uk[k] 
      lk[k] <- uk[k] <- NA
      
      
    }
  }
  
  return(list(uk=uk,
              lk=lk,
              ck=max(ck.unrestricted,cMin),
              ck.unrestricted=ck.unrestricted,
              alphaSpent=alphaSpent,
              betaSpent=betaSpent))
  
}
