#Notes:
#cannot reject H0 after recruitment is stopped due to futility boundary crossing at interim
#requires prediction of information at interim, but ck can be adjusted upon observing this info

#' @title Calculate boundaries for a group sequential design
#' @description Calculate boundaries for a group sequential design with delayed endpoints based on planned and/or observed information using an error spending approach. The function gives the boundaries for a non-binding futility rule using Method 2 as proposed by Hampson and Jennison.
#' @noRd
#' 
#' 
#' @param rho_alpha rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
#' @param rho_beta rho parameter of the rho-family spending functions (Kim-DeMets) for beta
#' @param alpha type I error
#' @param beta type II error
#' @param kMax max number of analyses (including final)
#' @param Info.max maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
#' @param InfoR.i Expected or observed (wherever possible) information rates at the interim and final analyses 1:Kmax
#' @param InfoR.d Expected or observed information rates at all potential decision analyses 1:(Kmax-1)
#' @param delta expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
#' @param abseps tolerance for precision when finding roots or computing integrals
#' @param direction greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
#' @param Trace Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
#' @param nWhileMax Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
#' @param toldiff Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
#' @param tolcoef Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
#' @param mycoefMax Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
#' @param mycoefL Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
#' @param myseed seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
#' @param sided one or two sided
#'  
#'
#' @examples
#'
#' Example to check that code matches
#' 
#' to reproduce bounds from CJ DSBS course slide 106
#'               
#' bCJ <- NonBindingHJ(rho_alpha=1.345,
#'            rho_beta=1.345,
#'            alpha=0.025,
#'            beta=0.1,
#'            Kmax=3,
#'            Info.max=12,
#'            InfoR.i=c(3.5,6.75,12)/12,
#'            InfoR.d=c(5.5,8.75)/12,
#'            delta=1,  
#'            abseps = 1e-06, 
#'            direction="smaller",
#'            sided=1
#'            )


#-------------
# No no longer need of this function, since binding option included in Method3 !!!!
#-----------

## * Method2_PC (code)
NonBindingHJ <- function(rho_alpha=2,        # rho parameter of the rho-family spending functions (Kim-DeMets) for alpha
                         rho_beta=2,           # rho parameter of the rho-family spending functions (Kim-DeMets) for beta
                         alpha=0.025,          # Type-I error (overall)
                         beta=0.2,             # Type-II error (overall)
                         Kmax,                 # number of planned analyses (including the final analysis)
                         Info.max=NULL,        # Info.max, i.e. maximum information needed for given beta (type II-error), delta (expected difference), alpha (type I-error) and Kmax. It can be given if it is known. Otherwise it is computed from the  values given for alpha, beta, delta and Kmax.
                         InfoR.i=NULL,         # Expected or observed (wherever possible) information rates at the interim and final analyses 1:Kmax
                         InfoR.d=NULL,         # Expected or observed information rates at all potential decision analyses 1:(Kmax-1)
                         delta=0,              # expected effect under the alternative (should be on the scale of the test statistc for which If and Info.max relate to one over the variance, e.g. delta=expected log(Hazard ratio))
                         abseps = 1e-06,       # tolerance for precision when finding roots or computing integrals
                         direction="smaller",  # greater is for Ho= theta > 0, "smaller" is for Ho= theta < 0 (note that in Jennison and Turnbull's book chapter (2013) they consider smaller)
                         Trace=FALSE,          # Used only if Info.max=NULL. Whether to print informations to follow the progression of the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax).
                         nWhileMax=30,         # Used only if Info.max=NULL. Maximum number of steps in the (root finding) algorithm to compute Info.max (from  alpha, beta, delta and Kmax)
                         toldiff= 1e-05,       # Used only if Info.max=NULL. Maximum tolerated difference between lower and upper bounds at anaylis Kmax (which souhld be zero), in the root finding algorithm, to find the value of Info.max
                         tolcoef= 1e-04,       # Used only if Info.max=NULL. Maximum tolerated difference before stopping the search (in the root finding algorithm), between two successive values for the multiplier coeficient 'coef' such that Info.max=coef*If  (some values for coef are given in Table 7.6 page 164 Jennison's book. The value of "If" (which stands for Information for Fixed design) corresponds to the information we need if Kmax=1)
                         mycoefMax= 1.2,       # Used only if Info.max=NULL. Upper limit of the interval of values in which we search for the multiplier coeficient 'coef' such that Info.max=coef*If (in the root finding algorithm).
                         mycoefL=1,            # Used only if Info.max=NULL. Lower limit of the interval (see mycoefMax)
                         myseed=2902,          # seed for producing reproducible results. Because we call functions which are based on Monte-Carlo compuation (pmvnorm)
                         sided=1              # one or two sided
                         ){
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
    cMin <- qnorm(1-alpha)
    mycoef <- NULL # initialize as needed for output
    lk <- rep(-Inf,Kmax)
    uk <- rep(Inf,Kmax) 
    ck <- rep(cMin,Kmax)
  
  if( (direction=="smaller" & delta<0) | (direction=="greater" & delta>0)){
    stop("The values given for arguments direction and delta are inconsistent.\n When direction=smaller, delta should be positive.\n When direction=greater, delta should be negative.")
  }
  
  if(direction=="greater"){
    delta <- -delta
  }else if(!direction%in%c("greater","smaller")){
      stop("direction should be either greater or smaller")
  }
  
  # initialize
  thealpha <- rep(0,Kmax)  # alpha spent up to step k
  thebeta <- rep(0,Kmax)   # beta spent up to step k
  IncAlpha <- rep(0,Kmax)  # alpha spent at step k
  IncBeta <- rep(0,Kmax)   # beta spent at step k     
 
  # compute variance-covariance matrix of vector (Z_1,...,Z_k)
  sigmaZk <- diag(1,Kmax)
  for(i in 1:Kmax){
    for(j in i:Kmax){
      sigmaZk[i,j] <- sqrt(InfoR.i[i]/InfoR.i[j])
      sigmaZk[j,i] <- sqrt(InfoR.i[i]/InfoR.i[j])
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
    xx <- NonBindingHJ(rho_alpha=rho_alpha,
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
                     direction=direction,
                     Trace=FALSE,
                     sided=sided)
    thediff <- abs(xx$boundaries[Kmax,"u.k"]-xx$boundaries[Kmax,"l.k"])
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
        xx <- NonBindingHJ(rho_alpha=rho_alpha,
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
                         direction=direction,
                         sided=sided)
        thediff <- abs(xx$boundaries[Kmax,"u.k"]-xx$boundaries[Kmax,"l.k"])
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
          ## browser()
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
    mycoef <- Info.max/If
  }
  ## {{{ compute vector of means under the alternative H1
  
  #compute information at each analysis
  Info.i <- InfoR.i*Info.max
  Info.d <- InfoR.d*Info.max
  
  # compute the mean of the multivariate normal distribution under the alternative H1
  thetheta <- delta*sqrt(Info.i)
  ## }}} 
  ## {{{ case k=1 
  IncAlpha[1] <- ErrorSpend(I=Info.i[1],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max)
  IncBeta[1] <-  ErrorSpend(I=Info.i[1],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)
  
  #efficacy boundary
  find.uk <- function(x){
    
    #information matrix for first interim and decision analysis
    sigmaZk2 <- matrix(NA,ncol=2,nrow=2)
    sigmaZk2[1,1] <- sigmaZk[1,1]
    sigmaZk2[2,2] <- 1
    sigmaZk2[1,2] <- sigmaZk2[2,1] <- sqrt(Info.i[1]/Info.d[1])
    
    #probability to conclude efficacy
    pmvnorm(lower = c(x,cMin),
            upper = c(Inf,Inf),
            mean=c(0,0),
            sigma= sigmaZk2,
            abseps = abseps) - IncAlpha[1]
  }
  
  uk[1] <- uniroot(find.uk,lower=-10,upper=10)$root  #dirty solution to use -10 and 10 for bounds
  
  #futility boundary
  find.lk <- function(x){
    
    #information matrix for first interim and decision analysis
    sigmaZk2 <- matrix(NA,ncol=2,nrow=2)
    sigmaZk2[1,1] <- sigmaZk[1,1]
    sigmaZk2[2,2] <- 1
    sigmaZk2[1,2] <- sigmaZk2[2,1] <- sqrt(Info.i[1]/Info.d[1])
    
    #probability to conclude futility
    pmvnorm(lower = c(uk[1],-Inf),
            upper = c(Inf,cMin),
            mean=c(thetheta[1],delta*sqrt(Info.d[1])),
            sigma= sigmaZk2,
            abseps = abseps) +
      pmvnorm(lower = c(-Inf,-Inf),
              upper = c(x,Inf),
              mean=c(thetheta[1],delta*sqrt(Info.d[1])),
              sigma= sigmaZk2,
              abseps = abseps) - IncBeta[1]
    
  }
  lk[1] <- uniroot(find.lk,lower=uk[1]-10,upper=uk[1])$root  #dirty solution to use -10 for lower bound
  
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
        thealpha[k] <- ErrorSpend(I=Info.i[k],rho=rho_alpha,beta_or_alpha=alpha,Info.max=Info.max) 
        IncAlpha[k] <- thealpha[k] - thealpha[(k-1)]   
        thebeta[k] <- ErrorSpend(I=Info.i[k],rho=rho_beta,beta_or_alpha=beta,Info.max=Info.max)  
        IncBeta[k] <- thebeta[k] - thebeta[(k-1)]   
        ## {{{ 
        ## {{{ u_k by solving what follows 
        uk[k] <- (uk[k-1] + lk[k-1])/2 # just to handle cases in which there is no root in what follows (when binding = TRUE )
        
        if(!k==Kmax){
          
          #information matrix for interim analyses and decision analysis
          sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
          sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
          sigmaZk2[k+1,k+1] <- 1
          sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.d[k])
          
          try(uk[k] <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf,k-1),x,cMin),
                                                   upper = c(uk[1:(k-1)],Inf,Inf),
                                                   mean=rep(0,k+1),
                                                   sigma= sigmaZk2,
                                                   abseps = abseps) - IncAlpha[k]},
                               lower = lk[k-1],
                               upper = uk[k-1],
                               tol = abseps)$root, silent = TRUE)
        } else {
          try(uk[k] <- uniroot(function(x){pmvnorm(lower = c(rep(-Inf,k-1),x),
                                                   upper = c(uk[1:(k-1)],Inf),
                                                   mean=rep(0,k),
                                                   sigma= sigmaZk[1:k,1:k],
                                                   abseps = abseps) - IncAlpha[k]},
                               lower = lk[k-1],
                               upper = uk[k-1],
                               tol = abseps)$root, silent = TRUE)
        }
        
        IsbkOK <- !(uk[k]==((uk[k-1] + lk[k-1])/2))
        
        if(!IsbkOK){warning(paste0("Could not compute uk[",k,"]"))}
        
        ## }}}
        ## {{{ a_k by solving what follows 
        if(IsbkOK){
           if(k!=Kmax){
             find.lkk <- function(x){
               
               #information matrix for interim analyses and decision analysis
               sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
               sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
               sigmaZk2[k+1,k+1] <- 1
               sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.d[k])
               
               #probability to conclude futility
               pmvnorm(lower = c(lk[1:(k-1)],uk[k],-Inf),
                       upper = c(uk[1:(k-1)],Inf,cMin),
                       mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                       sigma= sigmaZk2,
                       abseps = abseps) +
               pmvnorm(lower = c(lk[1:(k-1)],-Inf,-Inf),
                       upper = c(uk[1:(k-1)],x,Inf),
                       mean=c(thetheta[1:k],delta*sqrt(Info.d[k])),
                       sigma= sigmaZk2,
                       abseps = abseps) - IncBeta[k]
               
             }
             lk[k] <- try(uniroot(find.lkk,lower=uk[k]-10,upper=uk[k])$root)
              
             if(inherits(lk[k], "try-error")){
               lk[k] <- uk[k] # just to handle cases in which there is no root
               warning(paste0("try-error for calculation of lk[",k,"]"))
             }
           }
           if(k==Kmax){
             lk[k] <- uk[k] # just to handle cases in which there is no root
              
             try(lk[k] <- uniroot(function(x){pmvnorm(lower = c(lk[1:(k-1)],-Inf),
                                                      upper = c(uk[1:(k-1)],x),
                                                      mean=thetheta[1:k],
                                                      sigma= sigmaZk[1:k,1:k],
                                                      abseps = abseps) - IncBeta[k]},
                                                      lower = lk[k-1], 
                                                      upper = uk[k], 
                                                      tol = abseps)$root, silent = TRUE)
            if(inherits(lk[k],"try-error")){warning("try-error for calculation of lk[Kmax]")}
              
          }
            #----------------------------
            
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
  
  if(direction=="greater"){
    lk <- -lk
    uk <- -uk
    delta <- -delta
  }
  
  ck[Kmax] <- NA
  #ck[kMax] <- uk[kMax]
  #lk[kMax] <- uk[kMax] <- NA
  
  ## {{{ create  output
  d <- data.frame(l.k=lk,
                  u.k=uk,
                  c.k=ck,
                  Type.I.Error=thealpha,
                  Type.II.Error=thebeta,
                  Inc.Type.I=IncAlpha,
                  Inc.Type.II=IncBeta,
                  I.k=Info.i
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
              coef=mycoef,
              abseps=abseps,
              toldiff=toldiff,
              direction=direction,
              sided=sided,
              binding=FALSE,
              cMin=cMin)
  ## class(out) <- "delayedGSD"
  ## }}}
  .Random.seed <<- old # restore the current seed (before the call to the function)
  out
}
