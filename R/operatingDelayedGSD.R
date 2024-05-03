### operatingDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 29 2024 (09:41) 
## Version: 
## Last-Updated: maj  3 2024 (13:57) 
##           By: Brice Ozenne
##     Update #: 222
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * operatingDelayedGSD (documentation)
##' @title Operating Characteristics of GSD via simulation
##' @description Simulate trials with delayed endpoint and output the operating characteristics of a GSD.
##'
##' @param n.obs [integer vector, >0] total sample size in the trial.
##' If \code{NULL}, the sample size corresponds to the type II error of the first method will be used.
##' @param n.sim [integer vector, >0] number of simulations.
##' Can also be of length 2, in that case the second element indicates which replicate of the simulation is being run.
##' @param method [data.frame] method(s) to be used.
##' Should be a data.frame containing 3 columns: method (integer from 1 to 3), binding (logical), fixC (logical).
##' An additional column \code{overrule.futility} can be added to indicate what to do in presence of non-binding futility rules (default: never stop for futility under the null, always stop for futility under the alternative).
##' @param kMax [integer] max number of analyses (including final).
##' @param InfoR.i [numeric vector of length kMax] planned information rates at interim and final.
##' @param InfoR.d [numeric vector of length kMax] planned information rates at decision and final.
##' @param delta [numeric] planned treatment effect (i.e. difference in mean) on primary outcome.
##' 
##' @param args.GenData [list] arguments passed to the function GenData used to generate datasets. 
##' @param alpha [numeric, 0-1] type I error (one sided).
##' @param beta [numeric, 0-1] type II error (one sided).
##' @param rho_alpha [numeric, >0] rho parameter for alpha error spending function.
##' @param rho_beta [numeric, >0] rho parameter for beta error spending function.
##' @param PropForInterim [numeric vector of length Kmax-1] percentage of all subjects, once they had the change to have complete follow-up, triggering interim and final analyses.
##' By default equal to \code{InfoR.i}.
##' @param lag [numeric] time lag between stop of recruitment and decision to stop recruitment.
##' @param path [character] directory where to export the results.
##' If \code{NULL}, will be exported in the current directory.
##' If \code{NA}, will not be exported.
##' @param name [character] start of the filename.
##' @param export.tempo [logical] should results be exported at each iteration in a temporary file.
##' @param trace [logical] When \code{TRUE}, ouput is generate in the console as the execution of the function progresses.
##' @param seed [integer, >0] Random number generator (RNG) state used when starting resampling. The same seed is used for all methods.
##' @param cpus [integer, >0] number of child-processes for parallel evaluation.
##'
##' @examples
##'
##' df.method <- data.frame(method = 1:3, binding = TRUE, fixC = c(FALSE,FALSE,TRUE))
##' 
##' MyMissProb <- matrix(c(0.04807692, 0.05769231, 0.00961538, 0.88461538),  nrow = 2, ncol = 2,
##'                      dimnames = list(c("V1 missing", "V1 not missing"),c("V2 missing", "V2 not missing")))
##' 
##' args.GenData <- list(rand.block = c(1, 1, 0, 0),
##'                      allsd = c(2.5, 2.1, 2.4),
##'                      mean0 = c(10, 0, 0),
##'                      delta = c(0, 0.5, 1)*1,
##'                      ar = 30,
##'                      cor.01.1 = -0.15,
##'                      cor.ij.1 = 0.68,
##'                      cor.0j.1 = -0.27,
##'                      MissProb = MyMissProb,
##'                      DigitsOutcome = 2,
##'                      TimeFactor = 42,
##'                      DigitsTime = 0)
##' 
##' ##### 2 stages #### 
##' res2stage <- operatingDelayedGSD(n.sim = 10, 
##'              method = df.method,
##'              args.GenData = args.GenData,
##'              kMax = 2, InfoR.i = c(0.56,1), InfoR.d = c(0.65, 1), delta = 1,
##'              PropForInterim = 0.5, lag = 21,
##'              seed = 1:10)
##'
##' ##### 3 stages ####
##' res3stage <- operatingDelayedGSD(n.sim = 10, 
##'              method = df.method[1,,drop=FALSE],
##'              args.GenData = args.GenData,
##'              kMax = 3, InfoR.i = c(0.40,0.65,1), InfoR.d = c(0.50,0.75, 1), delta = 1,
##'              PropForInterim = c(0.35,0.6), lag = 21,
##'              seed = 1:10)
##'

## * operatingDelayedGSD (code)
##' @export
operatingDelayedGSD <- function(n.obs = NULL, n.sim, 
                                method, kMax, InfoR.i, InfoR.d, delta, 
                                args.GenData = NULL,
                                alpha = 0.025, beta = 0.1, rho_alpha = 2, rho_beta = 2, 
                                PropForInterim = InfoR.i, lag = 0, 
                                path = NULL, name = "simGSD", export.tempo = TRUE,
                                trace = TRUE, seed = NULL, cpus = 1){


    ## ** normalize user input
    ## *** n.obs
    if(!is.null(n.obs) && length(n.obs)!=1){
        stop("Argument \'n.obs\' should either be NULL or have length 1. \n")
    }

    ## *** n.sim
    if(length(n.sim)==1){
        n.run <- 1
        multirun <- FALSE
    }else if(length(n.sim)==2){
        n.run <- n.sim[2]
        n.sim <- n.sim[1]
        multirun <- TRUE
    }else{
        stop("Argument \'n.sim\' must have length 1 or 2. \n")
    }

    ## *** seed
    if(all(is.na(seed))){
        seed <- NULL
    }else{
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
    }
    if(!is.null(seed) && length(seed)>1 && length(seed)!=n.sim){
        stop("If argument \'seed\' has length greater than 1, it should have length the number of requested simulation: ",n.sim,".\n")
    }

    ## *** method
    if(!is.data.frame(method)){
        stop("Argment \'method\' should be a data.frame. \n")
    }
    valid.colmethod <- c("method","binding","fixC")
    if(any(names(method) %in% valid.colmethod == FALSE) || any(valid.colmethod %in% names(method) == FALSE)){
        stop("Argment \'method\' should be a data.frame with columns \"method\", \"binding\", and \"fixC\". \n")
    }
    if(any(method$method %in% 1:3 == FALSE)){
        stop("Column \"method\" of argument \'method\' should take values in 1,2,3. \n")
    }    
    if(any(method$binding %in% c(TRUE,FALSE) == FALSE)){
        stop("Column \"binding\" of argument \'method\' should be TRUE or FALSE. \n")
    }    
    if(any(method$fixC %in% c(TRUE,FALSE) == FALSE)){
        stop("Column \"fixC\" of argument \'method\' should be TRUE or FALSE. \n")
    }    
    if(any(method$method==3 & method$fixC == FALSE)){
        stop("Incorrect argument \'method\': method 3 should always have fixC equal to TRUE. \n")
    }
    
    ## *** args.GenData
    if(any(names(args.GenData) %in% c("n","seed"))){
        stop("Argument \'args.GenData\' should not contain \"n\" or \"seed\". \n")
    }

    ## *** kMax
    ## *** alpha
    ## *** beta
    ## *** rho_alpha
    ## *** rho_beta
    ## *** InfoR.i
    if(length(InfoR.i)!=kMax){
        stop("Length of argument \'InfoR.i\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax i.e. ",kMax,".\n",sep="")
    }
    if(InfoR.i[kMax]!=1){
        stop("Argument \'InfoR.i\' at kMax should be 1. \n")
    }
    if(any(diff(InfoR.i)<0)){
        stop("Argument \'InfoR.i\' should be increasing. \n")
    }

    ## *** InfoR.d
    if(length(InfoR.d)!=kMax){
        stop("Length of argument \'InfoR.d\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax i.e. ",kMax,".\n",sep="")
    }
    if(InfoR.d[kMax]!=1){
        stop("Argument \'InfoR.d\' at kMax should be 1. \n")
    }
    if(any(diff(InfoR.d)<0)){
        stop("Argument \'InfoR.d\' should be increasing. \n")
    }

    ## *** delta
    if(length(delta)!=1){
        stop("Argument \'delta\' should have length 1. \n")
    }

    ## *** PropForInterim
    if(length(PropForInterim)==kMax && PropForInterim[kMax]==1){
        PropForInterim <- PropForInterim[-kMax]
    }else if(length(PropForInterim)!=(kMax-1)){
        stop("Length of argument \'PropForInterim\' incompatible with argument \'kMax\'. \n",
             "Should be of length kMax i.e. ",kMax-1,".\n",sep="")
    }
    if(any(diff(PropForInterim)<0)){
        stop("Argument \'PropForInterim\' should be increasing. \n")
    }

    ## *** cpus
    if(length(cpus)!=1){
        stop("Argument \'cpus\' should have length 1.\n ")
    }else if(identical(cpus,"all")){
        cpus <- parallel::detectCores()
    }else if(!is.numeric(cpus) || cpus <=0 || cpus %% 1 != 0){
        stop("Argument \'cpus\' should be a positive integer or \'all\'.\n ")
    }else if(cpus>1 && cpus>parallel::detectCores()){
        stop("Only ",parallel::detectCores()," CPU cores are available. \n")
    }
    
    ## ** initialization

    ## *** generative model
    formals.GenData <- formals(GenData)
    if("N.fw" %in% names(args.GenData) == FALSE){
        args.GenData$N.fw <- eval(formals.GenData$N.fw)
    }
    N.fw <- args.GenData$N.fw
    for(iArgs in c("rand.block","allsd","mean0","delta","ar","cor.01.1","cor.ij.1","cor.0j.1","TimeFactor")){
        if(iArgs %in% names(args.GenData) == FALSE){
            args.GenData[[iArgs]] <- eval(formals.GenData[[iArgs]])
        }
    }

    ## *** boundaries
    n.method <- NROW(method)
    e.bound <- lapply(1:n.method, function(iM){
        CalcBoundaries(kMax = kMax,  
                       alpha = alpha, 
                       beta = beta,  
                       InfoR.i = InfoR.i,
                       InfoR.d = InfoR.d,  
                       rho_alpha = rho_alpha,  
                       rho_beta = rho_beta,  
                       method = method$method[iM],  
                       cNotBelowFixedc = method$fixC[iM],
                       bindingFutility = method$binding[iM],
                       delta = delta)
    })

    ## *** sample size
    if(is.null(n.obs)){
        sd_n.obs <- utils::tail(args.GenData$allsd,1)*sqrt(1-utils::tail(args.GenData$cor.0j.1,1)^2)
        nFix.sample <- 2*2*(sd_n.obs/delta)^2*(qnorm(1-beta)-qnorm(alpha))^2 ## fixed design (extra factor 2* because the formula is per arm)
        if(!is.null(args.GenData$MissProb)){
            nFix.sample <- nFix.sample/(1-(args.GenData$MissProb[1,1]+args.GenData$MissProb[2,1]))
        }
        n.obs <- ceiling(e.bound[[n.method]]$planned$InflationFactor * nFix.sample)
    }

    ## *** simulations
    vec.sim <- ((n.run-1)*n.sim + 1):(n.run*n.sim) ## indices of all iterations for this replicate
    
    ## *** seed
    if(!is.null(seed)){
        if(!is.null(get0(".Random.seed"))){ ## avoid error when .Random.seed do not exists, e.g. fresh R session with no call to RNG
            old <- .Random.seed # to save the current seed
            on.exit(.Random.seed <<- old) # restore the current seed (before the call to the function)
        }else{
            on.exit(rm(.Random.seed, envir=.GlobalEnv))
        }
        if(length(seed)==1 && n.sim*n.run != 1){
            set.seed(seed)
            allseeds <- sample.int(n = 1e5, size = n.sim*n.run, replace=FALSE) #x=1:(.Machine$integer.max) seems to be maximal possible
        }else{
            allseeds <- rep(NA, n.sim*n.run)
            allseeds[vec.sim] <- seed
            
        }        
    }

    ## *** interim
    n.interim <- kMax-1

    ## *** overrule futility
    if("overrule.futility" %in% names(method) == FALSE){
        ## Non binding: never stop for futility when simulating under the null and always stop for futility when simulating under the alternative
        ## (then the observed rejection rate should match the nominal type 1 or type 2 error)
        method$overrule.futility <- all(abs(args.GenData$delta)<1e-10) & (method$binding == FALSE)
    }else if(any(method$binding == TRUE & method$overrule.futility == TRUE)){
        stop("Incorrect argument \'method\': column overrule.futility cannot be TRUE for binding futility rules. \n")
    }
    
    ## ** Welcome message
    if(trace){
        if(multirun){
            cat("\t\tSimulation study for Group Sequential trial number ",n.run,". \n",
                "\t\t(",n.sim," iterations ",vec.sim[1],":",vec.sim[n.sim],", sample size: ",paste(n.obs, collapse = ", "),"). \n\n",sep="")
        }else{
            cat("\t\tSimulation study for Group Sequential trial. \n",
                "\t\t(",n.sim," iterations, sample size: ",paste(n.obs, collapse = ", "),"). \n\n",sep="")
        }
    }

    ## ** Loop
    RES <- NULL # initialize results to save

    if (cpus == 1) {
        for(iSim in 1:n.sim){ ## iSim <- 1
            if(!is.null(seed)){
                iSeed <- allseeds[vec.sim[iSim]]
            }else{
                iSeed <- NULL
            }

            if(trace){
                if(multirun){
                    cat("simulation ",iSim,"(",vec.sim[iSim],")/",n.sim,": seed ",iSeed,"\n",sep="")
                }else{
                    cat("simulation ",iSim,"/",n.sim,": seed ",iSeed,"\n",sep="")
                }
            }

            ## *** Generate data
            iData <- GenData(n = n.obs, 
                             N.fw = N.fw,
                             rand.block = args.GenData$rand.block,
                             allsd = args.GenData$allsd,
                             mean0 = args.GenData$mean0,
                             delta = args.GenData$delta, 
                             ar = args.GenData$ar,
                             cor.01.1 = args.GenData$cor.01.1,
                             cor.ij.1 = args.GenData$cor.ij.1,
                             cor.0j.1 = args.GenData$cor.0j.1,
                             seed = iSeed,
                             MissProb = args.GenData$MissProb,
                             DigitsOutcome = args.GenData$DigitsOutcome,
                             TimeFactor = args.GenData$TimeFactor,
                             DigitsTime = args.GenData$DigitsTime
                             )$d

            iLs.res <- lapply(1:n.method, function(iM){  ## iM <- 1
                iOut <- try(runDelayedGSD(iData, boundaries = e.bound[[iM]], N.fw = N.fw, PropForInterim = PropForInterim, lag = lag, overrule.futility = method$overrule.futility[iM]))
                if(inherits(iOut,"try-error")){
                    return(NULL)
                }else{
                    return(cbind(iOut, seed = iSeed))
                }
            })

            ## *** save results
            RES <- rbind(RES,do.call(rbind, iLs.res))
            if(!is.null(path) && export.tempo){
                save(RES,file=file.path(path,paste0(name,"(tempo)-",iter_sim,".rda")))
            }
        }
        
    }else{

        ## start cluster
        cl <- parallel::makeCluster(cpus)
        on.exit(parallel::stopCluster(cl))
        ## link to foreach
        doSNOW::registerDoSNOW(cl)

        ## other things to export
        ## parallel::clusterExport(cl, varlist = "seqSeed", envir = environment())

        ## progress bar
        if(trace>0){
            pb <- utils::txtProgressBar(max = n.sim, style = 3)          
            progress <- function(n){utils::setTxtProgressBar(pb, n)}
            opts <- list(progress = progress)
        }else{
            opts <- list()
        }
        

        iX <- NULL ## [:forCRANcheck:] foreach        
        iLs.res <- foreach::`%dopar%`(
                                foreach::foreach(iX=1:n.sim, .options.snow = opts, .packages = c("BB","emmeans","DelayedGSD","mvtnorm","nlme")), {
                                    if(!is.null(seed)){
                                        iSeed <- allseeds[vec.sim[iX]]
                                    }else{
                                        iSeed <- NULL
                                    }
                                    
                                    iData <- GenData(n = n.obs, 
                                                     N.fw = N.fw,
                                                     rand.block = args.GenData$rand.block,
                                                     allsd = args.GenData$allsd,
                                                     mean0 = args.GenData$mean0,
                                                     delta = args.GenData$delta, 
                                                     ar = args.GenData$ar,
                                                     cor.01.1 = args.GenData$cor.01.1,
                                                     cor.ij.1 = args.GenData$cor.ij.1,
                                                     cor.0j.1 = args.GenData$cor.0j.1,
                                                     seed = iSeed,
                                                     MissProb = args.GenData$MissProb,
                                                     DigitsOutcome = args.GenData$DigitsOutcome,
                                                     TimeFactor = args.GenData$TimeFactor,
                                                     DigitsTime = args.GenData$DigitsTime
                                                     )$d

                                    iLs.res <- lapply(1:n.method, function(iM){  ## iM <- 1
                                        iOut <- try(runDelayedGSD(iData, boundaries = e.bound[[iM]], N.fw = N.fw, PropForInterim = PropForInterim, lag = lag, overrule.futility = method$overrule.futility[iM]))
                                        if(inherits(iOut,"try-error")){
                                            return(NULL)
                                        }else{
                                            return(cbind(iOut, seed = iSeed))
                                        }
                                    })
                                    return(do.call(rbind, iLs.res))
                                })
        if(trace>0){close(pb)}        
        RES <- rbind(RES,do.call(rbind, iLs.res))
    }

    ## ** final export
    rownames(RES) <- NULL
    if(!is.null(path)){
        save(RES,file=file.path(path,paste0(name,"-",iter_sim,".rda")))
    }
    out <- list(call =  match.call(),
                args = c(list(n.obs = n.obs, n.sim = n.sim, kMax = kMax, method=method), args.GenData), 
                boundaries = e.bound,
                results = RES)
    class(out) <- append("operatingDelayedGSD",class(RES))
    return(out)
}   


##----------------------------------------------------------------------
### operatingDelayedGSD.R ends here
