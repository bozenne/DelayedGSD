#' @title Calculate boundaries for a group sequential design using Method 1
#' @description Calculate boundaries for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach.
#' 
#' 
#' 
#' @param rho_alpha rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
#' @param rho_beta rho parameter of the rho-family spending functions (Kim-DeMets) for beta
#' @param alpha type I error
#' @param beta type II error
#' @param kMax max number of analyses (including final)
#' @param Info.max maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
#' @param InfoR.i Expected or observed (wherever possible) information rates at the interim analyses 1:(Kmax-1)
#' @param InfoR.d Expected or observed information rates at all potential decision analyses including the final analysis 1:Kmax
#' @param delta expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
#' @param abseps tolerance for precision when finding roots or computing integrals
#' @param alternative a character string specifying the alternative hypothesis, \code{"greater"} or \code{"less"}.
#' H0 \eqn{\theta=0} vs H1 \eqn{theta<0} (\code{"less"}) or theta > 0 (\code{"greater"}).
#' Note that in Jennison and Turnbull's book chapter (2013) they consider greater.
#' @param binding whether we assume binding futility boundaries
#' @param Trace Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
#' @param nWhileMax Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
#' @param toldiff Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
#' @param tolcoef Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
#' @param mycoefMax Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
#' @param mycoefL Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
#' @param myseed seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
#' @param cMin minimun possible value c for the decision analysis, typically that for a fixed sample test (H & J page 10)
#' @param PowerCorrection whether or not to apply a correction to the type II error spending to reduce the extent to which Method 1 is overpowered   
#'  
#' @examples
#' ## see test-boundary.R                    


## * Method1 (code)
#' @export
Method1 <- function(rho_alpha=2,
                    rho_beta=2,
                    alpha=0.025,
                    beta=0.2,   
                    Kmax,       
                    Info.max=NULL, 
                    InfoR.i=NULL,
                    InfoR.d=NULL,
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
                    cMin=-Inf,
                    PowerCorrection=FALSE){
  
    require(mvtnorm)
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
    mycoef <- NULL # initialize as needed for output
    lk <- rep(-Inf,Kmax)
    uk <- rep(Inf,Kmax) 
    ck <- rep(NA,Kmax)
    ck.unrestricted <- rep(NA,Kmax)

    if(alternative=="greater" & delta<0){
        stop("The values given for arguments \'alternative\' and \'delta\' are inconsistent. \n",
             "When alternative=\"greater\", argument \'delta\' should be positive. \n")
    }else if(alternative=="less" & delta>0){
        stop("The values given for arguments alternative and delta are inconsistent. \n",
             "When alternative=\"less\", delta should be negative. \n")
    }
  
    if(alternative=="less"){
        delta <- -delta
    }else if(alternative != "greater"){
        stop("alternative should be either \"greater\" or \"less\".")
    }

    ## initialize
    thealpha <- rep(0,Kmax)  ## alpha spent up to step k
    thebeta <- rep(0,Kmax)   ## beta spent up to step k
    IncAlpha <- rep(0,Kmax)  ## alpha spent at step k
    IncBeta <- rep(0,Kmax)   ## beta spent at step k
    
    #information sequence relevant for alpha spending and covariance matrix
    InfoR <- c(InfoR.i,InfoR.d[Kmax])
    
    ## compute variance-covariance matrix of vector (Z_1,...,Z_k)
    sigmaZk <- diag(1,Kmax)
    for(i in 1:Kmax){
        for(j in i:Kmax){
            sigmaZk[i,j] <- sqrt(InfoR[i]/InfoR[j])
            sigmaZk[j,i] <- sqrt(InfoR[i]/InfoR[j])
        }
    }
                                        # compute If (see Jennison book page 87
    If <- (qnorm(1-alpha)+qnorm(1-beta))^2/delta^2
    if(Trace){
        cat("\n If computed as =",If,"\n")
    }
    ## }}}
    if(is.null(Info.max)){
                                        
        ## {{{ Compute Info.max from If and the other arguments (Recursive calls to the function)
        if(Trace){
            cat("\n We start the search of the value for coef=Info.max/If. \n \n")
        }        
        ## {{{ initialize key values to be updated in the following loop
        nwhile <- 0
        mycoefL0 <- mycoefL
        mycoefU <- mycoefMax
        mycoef <- mycoefU
        ## }}}
        ## {{{ Is the interval within which to search for coef large enough ?      
        if(Trace){
            cat("\n Check whether the interval within which to search for coef large enough. \n")
        }
        xx <- Method1(rho_alpha=rho_alpha,
                      rho_beta=rho_beta,
                      alpha=alpha,
                      beta=beta,
                      Kmax=Kmax,         
                      Info.max=If*mycoefU,
                      InfoR.i=InfoR.i,
                      InfoR.d=InfoR.d,
                      delta=delta,
                      abseps=abseps,
                      toldiff=toldiff,
                      alternative="greater",
                      binding=binding,
                      Trace=FALSE,
                      cMin=cMin,
                      PowerCorrection=PowerCorrection)
        thediff <- abs(xx$boundaries[Kmax,"uk"]-xx$boundaries[Kmax,"lk"])
        ## }}}
        if(thediff==0){
            ## {{{ if yes, we search coef
            mycoef <- (mycoefL + mycoefU)/2
            thediff <- 2*toldiff
            if(Trace){
                cat("\n we start the search within [",mycoefL,",",mycoefU,"] \n")
            }
            while(nwhile < nWhileMax & thediff>toldiff & abs(mycoefL-mycoefU)> tolcoef){
                nwhile <- nwhile + 1
                if(Trace){
                    cat("\n Step :",nwhile,"(out of max.", nWhileMax,")")
                }
                xx <- Method1(rho_alpha=rho_alpha,
                              rho_beta=rho_beta,
                              alpha=alpha,
                              beta=beta,
                              Kmax=Kmax,         
                              Info.max=If*mycoef,
                              InfoR.i=InfoR.i,
                              InfoR.d=InfoR.d,
                              delta=delta,
                              abseps=abseps,
                              toldiff=toldiff,
                              alternative="greater",
                              binding=binding,
                              cMin=cMin,
                              PowerCorrection=PowerCorrection)
                thediff <- abs(xx$boundaries[Kmax,"uk"]-xx$boundaries[Kmax,"lk"])
                if(thediff>toldiff){
                    if(Trace){
                        cat("\n Value coef=",mycoef,"is too small  \n")
                        cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
                        ## cat("\n b.K=",xx$boundaries[K,"b.k"]," and a.K=",xx$boundaries[Kmax,"a.k"]," \n")
                    }
                    mycoefL <- (mycoefL+mycoefU)/2
                    if(Trace){
                        cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
                    }
                    mycoef <- (mycoefL+mycoefU)/2
                }
                if(thediff==0){
                    if(Trace){
                        cat("\n Value coef=",mycoef,"is too large  \n")
                        cat("\n coef=",mycoef,"leads to b.K-a.K=",thediff, "(whereas tol=",toldiff,") \n")
                        ## cat("\n b.K=",xx$boundaries[Kmax,"b.k"]," and a.K=",xx$boundaries[Kmax,"a.k"]," \n")
                    }
                    mycoefU <- (mycoefL+mycoefU)/2
                    if(Trace){
                        cat("\n we update the interval : [",mycoefL,",",mycoefU,"] \n")
                    }
                    mycoef <- (mycoefL+mycoefU)/2
                    thediff <- 2*toldiff
                }
                if((thediff<=toldiff & thediff!=0) | abs(mycoefL-mycoefU)<= tolcoef ){
                    if(Trace){
                        cat("\n coef value FOUND : coef=",mycoef,"\n (leads to b.K-a.K=",thediff, " and tol.=",toldiff," and search interval length is=",abs(mycoefL-mycoefU),"and tol.=",tolcoef,")\n")
                    }
                    Info.max <- mycoef*If
                    if(Trace){
                        cat("\n Info.max computed as=",Info.max,"\n")
                    }
                    
                    ## print("Info.max created")
                }else{
                    if(nwhile==nWhileMax){
                        stop("Info.max could not be computed presicely enough : we need to allow for more iterations in the algorithm : you should probably call the function again with a larger value for nWhileMax.")
                    }
                }
            }
            ## }}}
        }else{
            ## {{{ if no, we stop and explain why
            stop("The interval [mycoefL,mycoefMax]= [",mycoefL0,",",mycoefMax,"] is too small. You should probably call the function again with a larger value for mycoefMax and/or a lower (value >=1) for mycoefL \n")
            ## }}}
        }        
    
    ## }}}
  }else{
    mycoef <- Info.max/If ## Inflation factor
  }
  ## {{{ compute vector of means under the alternative H1
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  Info <- InfoR*Info.max
  
  # compute the mean of the multivariate normal distribution under the alternative H1
  thetheta <- delta*sqrt(Info)
  ## }}} 
  ## {{{ case k=1 
  IncAlpha[1] <- ErrorSpend(I=Info[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
  IncBeta[1] <-  ErrorSpend(I=Info[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)

    uk[1] <- qnorm(p=1-IncAlpha[1],mean=0,sd=1)         # compute under the null (Ho)
    lk[1] <- qnorm(p=IncBeta[1],mean=thetheta[1],sd=1)  # compute under the alternative (H1)

    ck.unrestricted[1] <- calc_ck(uk=uk[1],
                                  lk=lk[1],
                                  Info.i=Info.i[1],
                                  Info.d=Info.d[1],
                                  Info.max=Info.max,
                                  ImaxAnticipated=FALSE,
                                  rho_alpha=rho_alpha,
                                  alpha=alpha,
                                  bindingFutility = binding)
    ck[1] <- max(ck.unrestricted[1], cMin)
                                        #------------
  
  thealpha[1] <- IncAlpha[1]   
  thebeta[1] <- IncBeta[1]
  
  lk <- pmin(lk,uk) # just in case of over-running
  
  ## if(Trace){
  ## cat("\n a.1 computed as",lk[1],"and b.1 as",uk[1]," \n")
  ## }
  ## }}}
  ## {{{ loop over k >=2
  
  
  if(Kmax>1){
    for(k in 2:Kmax){
      if(!lk[k-1]==uk[k-1]){
        ## {{{ if  over-running has not occurred yet
        thealpha[k] <- ErrorSpend(I=Info[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
        IncAlpha[k] <- thealpha[k] - thealpha[(k-1)]
        thebeta[k] <- ErrorSpend(I=Info[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)  
        if(PowerCorrection){
          #real type II error spent at analysis k-1
          IncBeta[k-1] <- TypeIIerrorSpent(lk=lk,
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
          thebeta[k-1] <- ifelse(k==2,IncBeta[k-1],thebeta[k-2]+IncBeta[k-1]) 
        } 
        IncBeta[k] <- thebeta[k] - thebeta[(k-1)]   
        
        
        ## {{{ 
        ## {{{ u_k by solving what follows 
        uk[k] <- (uk[k-1] + lk[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
        
        if(binding){
          TheLowerValues <- lk[1:(k-1)]
        }else{
          TheLowerValues <- rep(-Inf,k-1)
        }
        
        try(uk[k] <- uniroot(function(x){pmvnorm(lower = c(TheLowerValues,x),
                                                 upper = c(uk[1:(k-1)],Inf),
                                                 mean=rep(0,k),
                                                 sigma= sigmaZk[1:k,1:k],
                                                 abseps = abseps) - IncAlpha[k]},
                             lower = lk[k-1],
                             upper = uk[k-1],
                             tol = abseps)$root, silent = TRUE)
        IsbkOK <- !(uk[k]==((uk[k-1] + lk[k-1])/2))
        ## }}}
        ## {{{ a_k by solving what follows 
        if(IsbkOK){
          if(k!=Kmax){
            
            try(lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
                                                       upper = c(uk[1:(k-1)],x),
                                                       mean=thetheta[1:k],
                                                       sigma= sigmaZk[1:k,1:k],
                                                       abseps = abseps) - IncBeta[k]},
                                   lower = lk[k-1], 
                                   upper = uk[k], 
                                   tol = abseps)$root, silent = TRUE)
            
              if(!inherits(lk[k], "try-error")){
                  ck.unrestricted[k] <- calc_ck(uk=uk[1:k],
                                                lk=lk[1:k],
                                                Info.i=Info.i[1:k],
                                                Info.d=Info.d[k],
                                                Info.max=Info.max,
                                                ImaxAnticipated=FALSE,
                                                rho_alpha=rho_alpha,
                                                alpha=alpha,
                                                bindingFutility = binding)
                  ck[k] <- max(ck.unrestricted[k],cMin)
          } else {
              lk[k] <- uk[k] # just to handle cases in which there is no root
              if(inherits(lk[k],"try-error")){warning(paste0("try-error for calculation of lk[",k,"]"))}
          }
        }
          if(k==Kmax){
              lk[k] <- uk[k] # just to handle cases in which there is no root

              try(lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
                                                       upper = c(uk[1:(k-1)],x),
                                                       mean = thetheta[1:k],
                                                       sigma = sigmaZk[1:k,1:k],
                                                       abseps = abseps) - IncBeta[k]},
                                   lower = lk[k-1], 
                                   upper = uk[k-1],  ## put uk[k-1] instead of uk[k] (would often lead to no solution)
                                   tol = abseps)$root, silent = TRUE)
              if(inherits(lk[k],"try-error")){warning("try-error for calculation of lk[Kmax]")}
            
          }
        }else{
          lk[k] <- (uk[k-1] + lk[k-1])/2 # just to handle cases in which there is no root in what is above
        }
        ## }}}
        ## {{{ to deal with over-running (see chapter Jennison) 
        lk <- pmin(lk,uk)
        ## }}}
      }else{
        ## {{{ if  over-running has already occurred
        lk[k:Kmax] <- lk[k-1]
        uk[k:Kmax] <- uk[k-1]
        ## }}}
      }
    }
  }
  ## }}}
  
  ck.unrestricted[Kmax] <- uk[Kmax]
  ck[Kmax] <- max(uk[Kmax],cMin)
  
  if(alternative=="less"){
    lk <- -lk
    uk <- -uk
    delta <- -delta
    ck.unrestricted <- -ck.unrestricted
    ck <- -ck
  }
  
  #browser()
  
  ## {{{ create  output
  d <- data.frame(lk=lk,
                  uk=uk,
                  ck=ck,
                  ck.unrestricted=ck.unrestricted,
                  Type.I.Error=thealpha,
                  Type.II.Error=thebeta,
                  Inc.Type.I=IncAlpha,
                  Inc.Type.II=IncBeta,
                  Ik=Info
  )
  out <- list(boundaries=d,
              rho_alpha=rho_alpha,
              rho_beta=rho_beta,
              alpha=alpha,
              beta=beta,
              Kmax=Kmax,
              If=If,
              Info.max=Info.max,
              Info.i=Info.i,
              Info.d=Info.d,
              delta=delta,
              coef=mycoef, ## Inflation factor
              abseps=abseps,
              toldiff=toldiff,
              alternative=alternative,
              binding=binding,
              cMin=cMin,
              PowerCorrection=PowerCorrection)
  class(out) <- "delayedGSD"
  ## }}}
  return(out)
}
