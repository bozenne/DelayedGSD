#' @title Display Sample Size over Time
#' @description Display the available and pipeline sample size as a function of time
#'
#' @param d dataset
#' @param Delta.t ??
#' @param at ??
#' 
#' 
#' @export
PlotProgress <- function(d,Delta.t=0.500001,at=NA){
    n <- nrow(d) # final sample size
    allt <- unique(c(0,sort(as.matrix(d[,grep("t",names(d))]))))
    N.fw <- length(grep("t",names(d)))
    tmax <- max(d[,grep("t",names(d))],na.rm=T) + Delta.t # duration of trial (i.e. time of last collected data, eventually with additional delay Delta.t to do data management)
    allt <- c(allt,max(allt)+seq(from=0,to=1.5*Delta.t,by=min(diff(allt))))
    graphics::plot(x=c(0,tmax),y=c(0,n),col="white",xlab="Time since study start",ylab="Sample size",axes=FALSE)
    graphics::axis(1)
    graphics::axis(2,las=2)
    y0 <- sapply(allt,function(x) sum(d[,paste0("t",1)]<=x))
    graphics::lines(allt,y0,col=8,lwd=4)
    if(!is.na(at)){
        ny0 <- sum(d[,paste0("t",1)]<=at)
        graphics::text(x=at,y=ny0,pos=2,col=8,cex=2,labels=ny0)
    }
    for(i in 1:N.fw){
        y <- sapply(allt,function(x) sum(d[,paste0("t",i)]+Delta.t<=x))
        graphics::lines(allt,y,col=i,lwd=2)
        if(!is.na(at)){
            ny <- sum(d[,paste0("t",i)]+Delta.t<=at)
            graphics::text(x=at,y=ny,pos=2,col=i,cex=2,labels=ny)
        }        
    }
    for(i in 1:N.fw){
        yp <- sapply(allt,function(x) sum(d[,paste0("t",i)]+Delta.t>x & d[,"t1"]<=x))
        graphics::lines(allt,yp,col=i,lty=2,lwd=2)
        if(!is.na(at)){
            nyp <- sum(d[,paste0("t",i)]+Delta.t>at & d[,"t1"]<=at)
            graphics::text(x=at,y=nyp,pos=4,col=i,cex=2,labels=nyp)
        }
        
    }
    graphics::legend("left",title=paste0("Measurements data \n",ifelse(Delta.t==0,"collected","ready for analysis (plain) \n and in pipeline (dashed)")),fill=c(1:N.fw),legend=c(1:N.fw),bty="n")
    graphics::legend("topleft",col=8,legend="Subjects enrolled",bty="n",lwd=4)
    if(!is.na(at)){
        graphics::abline(v=at,lty=3)        
    }
    cbind(allt,y,yp)
}

#----------------------------------------------------------------------
### PlotProgress.R ends here
