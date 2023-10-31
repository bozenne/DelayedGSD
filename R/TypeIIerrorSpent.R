## * TypeIIerrorSpent
#Calculate the spent Type II error at a certain decision analysis which can be used to correct the type II error spending for Method 1
TypeIIerrorSpent <- function(lk,         #the futility boundaries (can be given only up to and including analysis k or can be complete boundaries)
                             uk,         #the efficacy boundaries (can be given only up to and including analysis k or can be complete boundaries)
                             ck,         #the decision boundaries (can be given only up to and including analysis k or can be complete boundaries)
                             Info.i,     #Expected or observed (wherever possible) information at the interim analyses 1:(Kmax-1) (can be given only up to analysis k)
                             Info.dk,    #Expected or observed (wherever possible) information at the decision analysis k
                             sigmaZk,    #The covariance matrix of the distribution of the test statistics (can be complete or only for interim analyses 1:k)
                             thetheta,   #The mean of the distribution of the test statistics under the alternative hypothesis
                             k,          #The analysis at which to calculate the spent Type II error
                             delta,      #The effect for which the study is powered
                             abseps){    #tolerance for precision when finding roots or computing integrals
  #browser()
  #construct covariance matrix for analysis k, including information at decision analysis k
  sigmaZk2 <- matrix(NA,ncol=k+1,nrow=k+1)
  sigmaZk2[1:k,1:k] <- sigmaZk[1:k,1:k]
  sigmaZk2[k+1,k+1] <- 1
  sigmaZk2[1:k,k+1] <- sigmaZk2[k+1,1:k] <- sqrt(Info.i[1:k]/Info.dk)
  
  #probability to stop for efficacy and conclude futility when there is an effect
  pmvnorm(lower = c(lk[0:(k-1)],uk[k],-Inf),
          upper = c(uk[0:(k-1)],Inf,ck[k]),
          mean=c(thetheta[0:k],delta*sqrt(Info.dk)),
          sigma= sigmaZk2,
          abseps = abseps) +
    #probability to stop for futility and conclude futility when there is an effect
    pmvnorm(lower = c(lk[0:(k-1)],-Inf,-Inf),
            upper = c(uk[0:(k-1)],lk[k],ck[k]),
            mean=c(thetheta[0:k],delta*sqrt(Info.dk)),
            sigma= sigmaZk2,
            abseps = abseps)
}
