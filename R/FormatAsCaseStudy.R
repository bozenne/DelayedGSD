#' @title Convert to L-format
#' @description Rename variable and relevel factor variables.
#' 
#' @param [data.frame] data output by \code{GenData}.
#' 
#' @export
FormatAsCase <- function(d){
    rm.col <- setdiff(names(d),paste0("missing",1:3))
    long <- stats::reshape(d[,rm.col,drop=FALSE],
                    direction='long',
                    varying=c('X2','X3'),
                    idvar='id',
                    v.names='Change')
    long <- long[order(long$id),]
    ## head(d)
    ## head(long)
    #----
    dd <- long
    names(dd) <- c("USUBJID","TRT01P","BASE","RANDDT","DT14","DT28","AVISIT","CHG")
    dd$TRT01P <- factor(dd$TRT01P,levels=c(0,1),labels=c("Control","Active"))
    dd$AVISIT <- factor(dd$AVISIT,levels=c(1,2),labels=c("Day 14","Day 28"))
    # time unit in days: outdated, now option in the main GenData function
    ## dd$RANDDT <- round(dd$RANDDT*14)
    ## dd$DT14 <- round(dd$DT14*14)
    ## dd$DT28 <- round(dd$DT28*14)
    # round with 2 digits: outdated, now option in the main GenData function
    ## dd$BASE <- round(dd$BASE,2)
    ## dd$CHG <- round(dd$CHG,2)
    ## head(dd)
    # Add FASFL
    IDFASFLN <- d[is.na(d$X2) & is.na(d$X3),"id"]
    dd$FASFL <- factor(dd$USUBJID %in% IDFASFLN,levels=c(TRUE,FALSE),labels=c("N","Y"))
    ## summary(dd)
    dd
}
