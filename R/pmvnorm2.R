### pmvnorm2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: okt 11 2023 (13:48) 
## Version: 
## Last-Updated: okt 11 2023 (13:50) 
##           By: Brice Ozenne
##     Update #: 3
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * pmvnorm2 (code)
pmvnorm2 <- function(lower, upper, mean, sigma, digits = 9){

    out <- mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean, sigma = sigma)

    if(is.na(out)){
        ## Handle the case where pmvnorm returns NaN, typically this is when integrating a domain with very little density.
        ## Sometimes just rounding the input help avoid NaN (which is strange)
        out <- mvtnorm::pmvnorm(lower = round(lower,digits), upper = round(upper,digits), mean = round(mean,digits), sigma = round(sigma,digits))
        ## Otherwise add-hoc criteria anticipating very little density
        if(is.na(out) && any(upper+5 < mean)){
            out <- 0
            attr(out,"msg") <- "pmvnorm returns NA - set to 0  (integration region way above the mean)"
        }else if(is.na(out) && any(lower - 5 > mean)){
            out <- 0
            attr(out,"msg") <- "pmvnorm returns NA - set to 0 (integration region way below the mean)"
        }else if(any(upper[1:(length(upper)-1)] + 2 < lower[2:length(lower)])){
            out <- 0
            attr(out,"msg") <- "pmvnorm returns NA - set to 0 (no overlap between consecutive integration regions)"
        }
    }

    return(out)
}


##----------------------------------------------------------------------
### pmvnorm2.R ends here
