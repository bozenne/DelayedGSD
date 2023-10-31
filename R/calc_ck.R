#' @title Calculate Critval c Using Method 1
#' @title calculate critval c based on the boundaries at interim.
#'
#' @param uk upper bounds for all analyses up to and including current stage k
#' @param lk lower bounds for all analyses up to and including current stage k
#' @param Info.i Information for all analyses up to and including current stage k
#' @param Info.d Observed information at decision analysis k
#' @param Info.max Maximum information
#' @param ImaxAnticipated set to TRUE if c should be calculated according to Eq 15 in HJ because study was stopped early due to Imax reached
#' @param rho_alpha value of rho used for type I error spending function
#' @param alpha one-sided alpha level to be used for the study
#' @param bindingFutility whether the futility stopping rule is binding
#'
#' @export
calc_ck <- function(uk,  
                    lk,  
                    Info.i,
                    Info.d,
                    Info.max=NULL,
                    ImaxAnticipated=FALSE,
                    rho_alpha=2, 
                    alpha=0.025, 
                    bindingFutility=TRUE){  
  
    requireNamespace("mvtnorm")
    requireNamespace("BB")
  
    ## message("the method assumes that positive effects are good")
    
    k <- length(uk)
    if(Info.d < Info.i[k]){
      stop("Decreasing information, cannot compute ck because covariance matrix not positive semidefinite")
    }
    
    ## If interim analysis k was skipped due to max info reached at interim or expected to be reached at decision, then calculate c according to Eq 15 from HJ
    if(ImaxAnticipated==TRUE){
        if(k==1){ ## we just have a single analysis of the trial where we spend all alpha
            return(qnorm(p=1-alpha,mean=0,sd=1))
        }
        
        Info.all <- c(Info.i[1:(k-1)],Info.d)
        sigmaZk <- diag(1,k)
        for(i in 1:(k)){
            for(j in i:(k)){
                sigmaZk[i,j] <- sqrt(Info.all[i]/Info.all[j])
                sigmaZk[j,i] <- sqrt(Info.all[i]/Info.all[j])
            }
        }
      
        if(bindingFutility){
            f <- function(x){
                y <- pmvnorm2(lower = c(lk[0:(k-1)],x),
                              upper = c(uk[0:(k-1)],Inf),
                              mean=rep(0,k),
                              sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.all[k-1], rho = rho_alpha,beta_or_alpha = alpha, Info.max = Info.max))
                as.numeric(y)
            }
        } else {
            lk2 <- rep(-Inf,length(lk))
        
            f <- function(x){
                y <- pmvnorm2(lower = c(lk2[0:(k-1)],x),
                              upper = c(uk[0:(k-1)],Inf),
                              mean=rep(0,k),
                              sigma= sigmaZk) - (alpha-ErrorSpend(I = Info.all[k-1], rho = rho_alpha,beta_or_alpha = alpha, Info.max = Info.max))                
                as.numeric(y)
            }
        }
        ## usual way of computing c based on Eq 14 HJ  
    } else {
      
      Info.all <- c(Info.i,Info.d)
      sigmaZk <- diag(1,k+1)
      for(i in 1:(k+1)){
        for(j in i:(k+1)){
          sigmaZk[i,j] <- sqrt(Info.all[i]/Info.all[j])
          sigmaZk[j,i] <- sqrt(Info.all[i]/Info.all[j])
        }
      }
      
        if(bindingFutility){
            f <- function(x){
                y <- pmvnorm2(lower = c(lk[0:(k-1)],uk[k],-Inf),
                              upper = c(uk[0:(k-1)],Inf,x),
                              mean = rep(0,k+1),
                              sigma = sigmaZk) - 
                    pmvnorm2(lower = c(lk[0:(k-1)],-Inf,x),
                             upper = c(uk[0:(k-1)],lk[k],Inf),
                             mean = rep(0,k+1),
                             sigma = sigmaZk)
                as.numeric(y)
            }
        } else {
            lk2 <- rep(-Inf,length(lk))
        
            f <- function(x){
                y <- pmvnorm2(lower = c(lk2[0:(k-1)],uk[k],-Inf),
                              upper = c(uk[0:(k-1)],Inf,x),
                              mean=rep(0,k+1),
                              sigma= sigmaZk) - 
                    pmvnorm2(lower = c(lk2[0:(k-1)],-Inf,x),
                             upper = c(uk[0:(k-1)],lk[k],Inf),
                             mean=rep(0,k+1),
                             sigma= sigmaZk)
                as.numeric(y)
            }
        }
    }

    lowerRoot <- lk[utils::tail(intersect(which(!is.infinite(lk)),1:k),1)]  ## last boundary among the k-1 already computed that is not infinite
    upperRoot <- uk[utils::tail(intersect(which(!is.infinite(uk)),1:k),1)]*1.1
    
    if(!is.null(Info.max)){
        if(ImaxAnticipated & Info.d > Info.max){
            lowerRoot <- -10
        }
    }

    if(is.nan(f(lowerRoot))){
        ## cannot evaluate the probability at lowerRoot as it is too small
        ccc <- BB::dfsane(f,
                          par = (upperRoot*1.1+lowerRoot)/2,
                          control = list(maxit = 1000, tol = 1e-7), quiet = TRUE)$par
    }else{        
        ccc <- try(stats::uniroot(f,
                                  lower = lowerRoot,
                                  upper = upperRoot)$root,
                   silent=TRUE)
        if(inherits(ccc,"try-error")){
            message("uniroot produced an error, switching to dfsane")
            ccc <- BB::dfsane(f,
                              par = (upperRoot*1.1+lowerRoot)/2,
                              control = list(maxit = 1000, tol = 1e-7), quiet = TRUE)$par
        }
    }
    return(ccc)
}
