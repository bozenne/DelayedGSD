## * updateMethod2 (documentation)

#' @title Update boundaries for a group sequential design using Method 2
#' @description Calculate boundaries for a group sequential design with delayed endpoints based observed information at interim/decision/final using an error spending approach.
#' @param rho_alpha rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
#' @param rho_beta rho parameter of the rho-family spending functions (Kim-DeMets) for beta
#' @param uk efficacy boundary from the previous interim analyses or planning
#' @param lk futility boundary from the previous interim analyses or planning
#' @param k [integer] current stage
#' @param type [character] type of analysis (\code{"interim"}, \code{"decision"}, \code{"final"}).
#' @param ImaxAnticipated [logical] Was the predicted information for this stage exceeding the maximum information?
#' @inheritParams Method2

## * updateMethod2 (code)
updateMethod2 <- function(rho_alpha=2,
                          rho_beta=2,
                          alpha=0.025,
                          alphaSpent = NULL,
                          beta=0.2,
                          betaSpent = NULL, 
                          Kmax,
                          Info.max=NULL, 
                          InfoR.i=NULL,  
                          InfoR.d=NULL,
                          uk = NULL,
                          lk = NULL,
                          k = NULL, type.k = NULL, ImaxAnticipated = FALSE,
                          delta=0,
                          abseps = 1e-06,
                          alternative="greater",
                          binding=TRUE,         
                          Trace=FALSE,          
                          nWhileMax=30,         
                          toldiff= 1e-05,       
                          tolcoef= 1e-04,       
                          mycoefMax= 1.2,       
                          mycoefL=1,            
                          myseed=2902,          
                          cMin=-Inf){

    ## {{{ set seed
    if(!is.null(myseed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }
        set.seed(myseed)
    }
    ## }}}
    ## {{{ preliminaries
    ## *** initialize boundaries
    ## lk <- rep(-Inf,Kmax)
    ## uk <- rep(Inf,Kmax)    
    ## initialize

    ## ## *** identify interim analyses that have been skipped due to decreasing information
    ##         index.DI <- which(object$conclusion["reason.interim",]=="decreasing information")
    ##         index.nDI <- setdiff(1:kMax,index.DI)
    ##         if(length(index.DI)>0){
    ##             ## When the interim analysis was skiped because of decreasing information
    ##             ## remove interim and decision
    ##             Info.i <- Info.i[index.nDI]
    ##             Info.d <- Info.d[setdiff(index.nDI,kMax)]
    ##             kMax <- kMax - length(index.DI) ## need to be after the update of InfoR.d as the previous line dependens on kMax
    ##         }



    ##         kMax <- kMax + length(index.DI)

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
                                        # compute If (see Jennison book page 87
    ## }}}
    ## {{{ compute vector of means under the alternative H1
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  Info <- InfoR*Info.max
  
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
                                        #thea[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)

                                        #------------
            find.lk <- function(x){
                                        #calculate c corresponding to lk
                ck <- max(calc_ck(uk=uk[1],
                                  lk=x,
                                  Info.i=Info.i[1],
                                  Info.d=Info.d[1],
                                  Info.max=Info.max,
                                  ImaxAnticipated=FALSE, ## Should not be here at interim with I>Imax because this case is handled in the parent function (updateBoundary)
                                  rho_alpha=rho_alpha,
                                  alpha=alpha,
                                  bindingFutility = binding), cMin)
    
                ## information matrix for first interim and decision analysis
                sigmaZk2 <- matrix(NA,ncol=2,nrow=2)
                sigmaZk2[1,1] <- sigmaZk[1,1]
                sigmaZk2[2,2] <- 1
                sigmaZk2[1,2] <- sigmaZk2[2,1] <- sqrt(Info.i[1]/Info.d[1])
    
                ## probability to conclude futility
                return(pmvnorm(lower = c(uk[1],-Inf),
                               upper = c(Inf,ck),
                               mean=c(thetheta[1],delta*sqrt(Info.d[1])),
                               sigma= sigmaZk2,
                               abseps = abseps) +
                       pmvnorm(lower = c(-Inf,-Inf),
                               upper = c(x,ck),
                               mean=c(thetheta[1],delta*sqrt(Info.d[1])),
                               sigma= sigmaZk2,
                               abseps = abseps) - betaSpent[1])
      
            }
            lk[1] <- uniroot(find.lk,lower=uk[1]-10,upper=uk[1])$root  #dirty solution to use -10 for lower bound
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
            find.lkk <- function(x){
                ## calculate c corresponding to lk
                ck <- max(calc_ck(uk=uk[1:k],
                              lk=c(lk[1:(k-1)],x),
                              Info.i=Info.i[1:k],
                              Info.d=Info.d[k],
                              Info.max=Info.max,
                              ImaxAnticipated=FALSE, ## Should not be here at interim with I>Imax because this case is handled in the parent function (updateBoundary)
                              rho_alpha=rho_alpha,
                              alpha=alpha,
                              bindingFutility = binding), cMin)
                
                ## information matrix for first interim and decision analysis
                sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
                sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
                sigmaZk2[k+1,k+1] <- 1
                sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.d[k])
                
                ## probability to conclude futility
                return(pmvnorm(lower = c(lk[1:(k-1)],uk[k],-Inf),
                               upper = c(uk[1:(k-1)],Inf,ck),
                               mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                               sigma= sigmaZk2,
                               abseps = abseps) +
                       pmvnorm(lower = c(lk[1:(k-1)],-Inf,-Inf),
                               upper = c(uk[1:(k-1)],x,ck),
                               mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                               sigma= sigmaZk2,
                               abseps = abseps) - betaSpentInc[k])
              
            }

            lk[k] <- uniroot(find.lkk,lower=uk[k]-10,upper=uk[k])$root
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
            betaSpent[k] <- beta
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
