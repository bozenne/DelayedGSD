#' @title Generate longitudnal (full) data
#' @description yy
#' 
#' @param n sample size
#' @param N.fw  number of follow-up measurements (equally spaced, after baseline)
#' @param rand.block for block randomization
#' @param allsd vector of sd of primary outcome at baseline and end of follow-up (main outcome)
#' @param mean0 mean outcome at each visit in control group
#' @param delta treatment effect (i.e. difference in mean) on primary outcome at each visit
#' @param ar accrual rate (average, unit is per time between the equally spaced visits)
#' @param cor.01.1 correlation between outcome at baseline and at first visit  (main outcome)
#' @param cor.ij.1 correlation between outcome at two consecutive follow-up measurements  (main outcome)
#' @param cor.0j.1 correlation between outcome at baseline and at any visit after the first visit  (main outcome)
#' @param seed integer for the random seed generator state.
#' @param MissProb Missingness probability. Should be an N.fw-dimensional array with in each dimension the proportion of missing and non missing (in that order).
#' For instance for two follow-up, should be a matrix containing the probability of missing both (1,1),
#' missing only the first timepoint (1,2)
#' missing only the second timepoint (2,1)
#' or have complete data (2,2).
#' @param DigitsOutcome Number of digits to round the outcome values (NULL means no rounding)
#' @param TimeFactor Multiply the times by a factor (e.g. 14 if time between two follow-up visit should be approx 14 days)
#' @param DigitsTime Number of digits to round the times (NULL means no rounding)
#' 
#' @details zz
#' @return ff
#' @author Paul Blanche
#' 
#' @examples
#' x <- GenData()
#' head(x$d,n=20)
#' 
#' @export
GenData <- function(n=52, 
                    N.fw=2,
                    rand.block=c(1,1,0,0),
                    allsd=rep(3,N.fw+1),
                    mean0=rep(0,N.fw+1),
                    delta=rep(0,N.fw+1),
                    ar=0.86*2,
                    cor.01.1=-0.15,
                    cor.ij.1=0.68,
                    cor.0j.1=-0.27,
                    seed=24082020,
                    MissProb=NULL,
                    DigitsOutcome=NULL,
                    TimeFactor=1,
                    DigitsTime=NULL
                    ){
    #--
    requireNamespace("mvtnorm")
    #--
    set.seed(seed)   
    d <- data.frame(id=1:n, # id
                    Z=blockrand(seed,rand.block,n))
    ## Z=c(rep(1,floor(n/2)),rep(0,n-floor(n/2)))[sample(1:n)])     # treatment
    NV <- N.fw + 1 # number of measurements (follow-up + baseline)
    for(k in 1:NV){
        d[[paste0("X",k)]] <- rep(NA,n) # initialize measurement values at each follow-up
        ## d[[paste0("t",k)]] <- rep(NA,n) # initialize time of each follow-up measurement
    }
    # correlation matrix of outcome at baseline and each follow-up time (can be made an independent subfunction)
    cor.matrix <- matrix(NA,ncol=NV,nrow=NV)
    diag(cor.matrix) <- 1
    cor.matrix[1,2] <- cor.matrix[2,1] <- cor.01.1
    if(NV>2){
        for(k in 3:NV){
            for(l in 1:(k-1)){
                cor.matrix[k,l] <- cor.matrix[l,k] <- cor.ij.1
            }
        }
        cor.matrix[1,3:NV] <- cor.matrix[3:NV,1] <- cor.0j.1
    }
                                        # Combine cor matrix and var vector to var-cov matrix
    vcovmat <- diag(allsd)%*%cor.matrix%*%diag(allsd)
    # Generate outcome at baseline and each follow-up visit
    allmeans1 <- delta + mean0  # mean in treatment group
    d[d$Z==0,paste0("X",1:NV)] <- mvtnorm::rmvnorm(n=sum(d$Z==0),mean=mean0,sigma=vcovmat)
    d[d$Z==1,paste0("X",1:NV)] <- mvtnorm::rmvnorm(n=sum(d$Z==1),mean=allmeans1,sigma=vcovmat)    
    # Generate time at inclusion (time since study start) and deduce time at each follow-up

    d[1,"t1"] <- 0
    for(k in 2:NV){            
        d[1,paste0("t",k)] <- d[1,"t1"] + k-1 # initialize time of each follow-up measurement
    }
    
    for(i in 2:n){
        d[i,"t1"] <- d[i-1,"t1"] +    stats::runif(1,min=0,max= 2/ar)
        for(k in 2:NV){            
            d[i,paste0("t",k)] <- d[i,"t1"] + k-1 # initialize time of each follow-up measurement
        }
    }

                                        # Missing values
    if(!is.null(MissProb)){        
        if(any(is.na(MissProb))){
            stop("Elements in argument \'MissProb\' should not be NAs. \n")
        }
        
        if(abs(sum(MissProb)-1)>1/(10*n)){
            stop("Elements in argument \'MissProb\' should sum up to one. \n")
        }
        if(length(dim(MissProb))!=N.fw){
            stop("Argument \'MissProb\' must have as many dimensions as number of follow-up (here ",N.fw,"). \n")
        }
        if(any(dim(MissProb)!=2)){
            stop("Each dimension of \'MissProb\' must have two elements: \"missing\" and \"non-missing\". \n")
        }
        dimnames.MissProb <- dimnames(MissProb)
        
        FreqMiss <- round(MissProb*n)
        FreqMiss[length(MissProb)] <- 0 ## last case is no missing value
        HowManyMiss <- sum(as.double(FreqMiss))
        whichMiss <- sample(x=1:n,size=HowManyMiss,replace=FALSE)
        
        typeMiss <- which(FreqMiss>0, arr.ind = TRUE)
        n.typeMiss <- FreqMiss[FreqMiss>0]
        nCum.typeMiss <- cumsum(n.typeMiss)
        
        if(length(whichMiss)>0){
            for(iType in 1:NROW(typeMiss)){ ## iType <- 3
                iX <- paste0("X",which(typeMiss[iType,]==1)+1)
                d[whichMiss[c(1,nCum.typeMiss+1)[iType]:nCum.typeMiss[iType]],iX] <- NA
            }
        }
    }


    ## round the outcome values
    if(!is.null(DigitsOutcome)){
        d[,paste0("X",1:NV)] <- round(d[,paste0("X",1:NV)],DigitsOutcome)
    }
    if(TimeFactor!=1){
        d[,paste0("t",1:NV)] <- d[,paste0("t",1:NV)]*TimeFactor
    }
    if(!is.null(DigitsTime)){
        d[,paste0("t",1:NV)] <- round(d[,paste0("t",1:NV)],DigitsTime)
    }   
    
    ## head(d)
    list(input=list(
             n=n, 
             N.fw=N.fw,
             rand.block=rand.block,
             allsd=allsd,
             mean0=mean0,
             delta=delta,
             ar=ar,
             cor.01.1=cor.01.1,
             cor.ij.1=cor.ij.1,
             seed=seed,
             MissProb=MissProb
         ),
         d=d)
}
