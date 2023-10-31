## * plot.delayedGSD (documentation)
#' @title Display boundaries over information
#'
#' @param x object from CalcBoundaries
#' @param y not used. For compatibility with the generic method.
#' @param type type of boundaries to plot (Z-statistic (Z), p-value (P), effect (E))
#' @param Itype information scale on x-axis (rate or absolute (abs))
#' @param main to specify a plot title
#' @param xlim range of the x axis
#' @param ylim range of the y axis
#' @param planned [logical] should the planned boundaries be displayed along with the boundaries that have been computed based on the data.
#' Can also be \code{"only"} to only display planned boundaries.
#' @param predicted [logical] Should the predicted information/boundaries at decision based on interim data be output (when relevant).
#' @param legend should a caption be added?
#' @param legend.x location of the legend on the plot.
#' @param legend.ncol number of column in the caption
#' @param legend.cex character expansion factor for the text in the legend.
#' @param add.arrow [logical] should arrows indicating the ordering of the space be added?
#' @param sep.arrow [numeric vector of length 2] space between the critical point and the arrow (below and above)
#' @param size.arrow [numeric, >0] width of the arrow.
#' @param lwd.arrow [numeric, >0] thickness of the arrow.
#' @param cex.bound [numeric, >0] size of the dots for the boundaries.
#' @param cex.estimate [numeric vector of size 2] size of the dots for the estimated statistic/p-value/effect and its label.
#' @param ... not used. For compatibility with the generic method.
#'
#' @examples
#' bnds <- CalcBoundaries(kMax=3, InfoR.i = c(0.5,0.75,1), InfoR.d = c(0.55,0.8))
#' plot(bnds, type = "Z")
#' plot(bnds, type = "Z", add.arrow = TRUE)
#' plot(bnds, type = "Z", add.arrow = TRUE, lwd.arrow=2,ylim=c(0.35,2.6), xlim=c(0.5,1.1),legend=FALSE)
#' plot(bnds, type = "E")
#' plot(bnds, type = "P")

## * plot.delayedGSD (code)
#' @export
plot.delayedGSD <- function(x,
                            y,
                            type="Z",   
                            Itype="rate",
                            planned=NULL,
                            predicted=TRUE,
                            main = NULL,   
                            xlim = NULL, 
                            ylim = NULL, 
                            legend=TRUE,  
                            legend.x=NULL,
                            legend.ncol=1,
                            legend.cex=1,
                            add.arrow=FALSE,
                            sep.arrow=c(0.2,0.1),
                            size.arrow=0.1,
                            lwd.arrow=1,
                            cex.bound=1.2,
                            cex.estimate=1.2,
                            ...){

                                        # {{{ Check user input
    if(!missing(y)){
        stop("Arugment \'y\' not used. \n")
    }
    dots <- list(...)
    if(length(dots)>0){
        stop("Unknown argument(s) \'",paste(names(dots),collapse="\' \'"),"\'. \n")
    }
    type <- match.arg(type, c("Z","E","P"))
    Itype <- match.arg(Itype, c("rate","abs"))
                                        # }}}
    
                                        # {{{ Preliminaries

    alpha <- x$alpha
    kMax <- x$kMax
    k <- x$stage$k
    kmain <- k
    type.k <- x$stage$type

    specdec <- function(x,k){ format(round(x,k),nsmall=k)}
    if(is.null(planned)){
        planned <- unlist(list(TRUE,"only")[(type.k=="planning")+1])
    }else if(planned %in% c(TRUE,FALSE) == FALSE && planned != "only"){
        stop("Argument \'planned\' should be TRUE, FALSE, or \"only\". \n")
    }
    
    
    ## extract from object
    if(identical(planned,"only")){
        outInfo <- coef(x, type = "information", planned = TRUE)
        outBound <- coef(x, type = "boundary", planned = TRUE)
        delta <- NULL
        index.skip <- NULL
    }else{
        outInfo <- coef(x, type = "information", planned = FALSE, predicted = predicted)
        outBound <- coef(x, type = "boundary", planned = FALSE, predicted = predicted)
        outBound <- outBound[,c("stage","Fbound","Ebound","Cbound")]
        outDelta <- stats::confint(x, k = "all", method = "ML")

        index.skip <- which(x$conclusion["reason.interim",]=="decreasing information")
        if(length(index.skip)>0){
            outDelta <- outDelta[outDelta$stage %in% index.skip == FALSE,,drop=FALSE]
            kmain <- k
            k <- k - length(index.skip)
        }
        delta <- switch(EXPR = type,
                        "Z" = outDelta$statistic,
                        "E" = outDelta$estimate,
                        "P" = outDelta$p.value)
        if(planned){
            outInfo2 <- coef(x, type = "information", planned = TRUE)
            outInfo[is.na(outInfo)] <- outInfo2[is.na(outInfo)]
            outBound2 <- coef(x, type = "boundary", planned = TRUE)
            outBound[is.na(outBound)] <- outBound2[is.na(outBound)]
        }
    }

    Info.i <- outInfo$Interim
    Info.d <- outInfo$Decision
    Info.max <- attr(outInfo,"Info.max")
   
    uk <- outBound[,"Ebound"]
    lk <- outBound[,"Fbound"]
    ck <- outBound[,"Cbound"]

    
    if(is.null(main)){
        if(identical(planned,"only")){
            main <- "Planning"
        }else if(type.k=="final"){
            main <- paste0("Final analysis (stage ",kmain,")")
        }else if(type.k=="decision"){
            main <- paste0("Decision analysis at stage ",kmain)
        }else if(type.k=="interim"){
            main <- paste0("Interim analysis at stage ",kmain)
        }
      
      if(length(index.skip)>0){
        main <- paste0(main, ", analysis ",paste(index.skip,collapse=",")," skipped")
      }
    }

    ## prepare
    if(Itype=="rate"){
        Ival <- Info.i/Info.max
        Idval <- Info.d/Info.max
        xlb <- "Information fraction"
    } else if(Itype=="abs"){
        Ival <- Info.i
        Idval <- Info.d
        xlb <- "Information"
    } 
    
    if(type=="P"){
        mydigits <- 3
    }else{
        mydigits <- 2
    }
                                        # }}}
                                        # {{{ Values to plot
                                        # default (if type="Z")

    ## NOTE: subset by [!is.na(uk)] to handle the case where we end the study early, e.g. for efficacy
    ## in that case bounds after the decision analysis are set to NA and should be ignored
    xu <- Ival[!is.na(uk) & !is.infinite(uk)]
    xd <- Idval[!is.na(ck) & !is.infinite(uk)]
    yu <- uk[!is.na(uk) & !is.infinite(uk)]
    yl <- lk[!is.na(lk) & !is.infinite(lk)]   
    yc <- ck[!is.na(ck) & !is.infinite(uk)]
    yh <- stats::qnorm(1-alpha)
    ylab <- "Stopping boundary (Z-statistic)"
 
    if(type=="P"){
        if(is.null(legend.x)){
            legend.x <- "topleft"
        }
        ylab <- "Stopping boundary (P-value)"
        yu <- 1-stats::pnorm(yu)
        yl <- 1-stats::pnorm(yl)
        yc <- 1-stats::pnorm(yc)
        yh <- alpha
    }else{
        if(is.null(legend.x)){
            legend.x <- "bottomright"
        }
    
    }
    
    if(type=="E"){
        ylab <- "Stopping boundary (Effect estimate)"
        yu <- yu/sqrt(Info.i)[!is.na(uk) & !is.infinite(uk)]
        yl <- yl/sqrt(Info.i)[!is.na(lk) & !is.infinite(lk)]
        yc <- yc/sqrt(Info.d)[!is.na(ck) & !is.infinite(uk)]
        yh <- stats::qnorm(1-alpha)/sqrt(Info.d[length(Info.d)])
    }
    if(is.null(ylim)){
        ylim <- switch(EXPR = type,
                       "P" = c(0,1),
                       "E" = c(min(c(yl[!is.infinite(yl)],delta),na.rm=TRUE)-0.5,max(c(yu[!is.infinite(yu)],delta),na.rm=TRUE)+0.5),
                       "Z" = c(min(c(-1,yl[!is.infinite(yl)],delta),na.rm=TRUE),max(c(3,yu[!is.infinite(yu)],delta),na.rm=TRUE))
                       )
    }
    if(is.null(xlim)){
        if(Itype=="rate"){
            xlim <- c(0,max(1,xd)*1.1)
        }else{
            xlim <- c(0,1.1*max(xd,na.rm=TRUE))
        }
    }
                                        # }}}
                                        # {{{ Plot
    graphics::plot(c(xu,xd[length(xd)]),c(yu,yc[length(yc)]),type="l",lty=2,lwd=2,ylim=ylim,col="green3",xlab=xlb,ylab=ylab,axes=FALSE,
                   xlim=xlim,main=main)
    graphics::lines(c(xu,xd[length(xd)]),c(yl,yc[length(yc)]),col="red",lwd=2,lty=2)
    graphics::points(xu,yl,col="red",pch=21,bg="red",cex=cex.bound)
    graphics::points(xu,yu,col="green3",pch=21,bg="green3",cex=cex.bound)
    graphics::points(xd,yc,col="black",pch=19,cex=cex.bound)

    if(!is.null(delta)){
        xdelta <- xu[1:k]
      
        if(type.k=="decision"){
            xdelta <- c(xdelta,xd[k])
        } else if (type.k=="final") {
          xdelta <- c(xu,xd[k])
        }
        ## lines(c(0,xdelta),c(0,delta),col="purple",lwd=2,lty=3)
        graphics::points(xdelta,delta,col="purple",pch=22,bg="purple",cex=cex.estimate[1])
        if(length(cex.estimate)>1){

            if(length(index.skip)>0){
                k.labels <- (1:x$stage$k)[-index.skip]
            }else{
                k.labels <- 1:k
            }

            if(type.k=="decision"){
                graphics::text(x=xdelta,y=delta,labels=c(k.labels,"D"),col="white", cex = cex.estimate[2])
            }else{
                graphics::text(x=xdelta,y=delta,labels=k.labels,col="white", cex = cex.estimate[2])
            }
        }
        graphics::text(x=xdelta,y=delta,labels=specdec(delta,k=mydigits),col="purple", pos = 4)
        
    }
    if(k>0){
        x.decision <- coef(x, type = "decision")
        if(!is.na(x.decision["comment",k]) && x.decision["comment",k]!="Imax reached"){
            graphics::text(x=xd,y=yc,labels=specdec(yc,k=mydigits),col="black",pos=2)
        }
    }
    
    if(planned == "only"){
        pos.l <- 1
        pos.u <- 3
    }else if(length(yl)==kMax){
        pos.l <- c(rep(2,k),rep(1,kMax-k))
        pos.u <- c(rep(2,k),rep(3,kMax-k))
    }else{
        pos.l <- 2
        pos.u <- 2
    }
    graphics::text(x=xu,y=yl,labels=specdec(yl,k=mydigits),col="red",pos=pos.l)
    graphics::text(x=xu,y=yu,labels=specdec(yu,k=mydigits),col="green3",pos=pos.u)    
    graphics::text(x=xd,y=yc,labels=specdec(yc,k=mydigits),col="black",pos=pos.l)
                                        #---   
    graphics::abline(h=yh,lty=2,col="grey")
    graphics::text(x=0,y=yh,labels=specdec(yh,k=mydigits),col="grey",pos=1)        
                                        # axes
    graphics::axis(1,at=c(0,xu,xd),labels=format(c(0,xu,xd),digits=2))
    graphics::axis(2,at=c(yl,yu,ylim),
                   labels=format(c(yl,yu,ylim),digits=2),las=2)
                                        # color areas
    colvectgreen <- grDevices::col2rgb("green3")/255
    colvectred <- grDevices::col2rgb("red")/255
    myrgbcolgreen <- grDevices::rgb(red=colvectgreen[1],
                                    green=colvectgreen[2],
                                    blue=colvectgreen[3],
                                    alpha=0.25)
    myrgbcolred <- grDevices::rgb(red=colvectred[1],
                                  green=colvectred[2],
                                  blue=colvectred[3],
                                  alpha=0.25)

    graphics::polygon(x=c(0,c(xu,xd[length(xd)]),rev(c(xu,xd[length(xd)])),0),
                      y=c(yl[1],c(yl,yc[length(yc)]),rep(ifelse(type=="P",max(ylim),min(ylim)),length(yu)+2)),
                      col=myrgbcolred,border=NA)
    graphics::polygon(x=c(0,c(xu,xd[length(xd)]),rev(c(xu,xd[length(xd)])),0),
                      y=c(yu[1],c(yu,yc[length(yc)]),rep(ifelse(type!="P",max(ylim),min(ylim)),length(yu)+2)),
                      col=myrgbcolgreen,border=NA)
                                        # legend   
    if(legend){
        if(length(legend.cex)==1){
            legend.cex <- rep(legend.cex,2)
        }
        
        text.legend <- c("Stopping bound for efficacy",
                         "Stopping bound for futility",
                         "Critical value fixed design",
                         "Decision boundary")
        lty.legend  <- c(2,2,2,NA)
        pch.legend  <- c(21,21,NA,19)
        col.legend  <- c("green3","red","grey","black")
        pt.bg.legend <- c("green3","red",NA,NA)
        if(!is.null(delta)){
            text.legend <- c(text.legend,"Estimate")
            lty.legend <- c(lty.legend,NA)
            pch.legend <- c(pch.legend,22)
            col.legend <- c(col.legend,"purple")
            pt.bg.legend <- c(pt.bg.legend,"purple")
        }
        
        graphics::legend(x = legend.x,
                         legend = text.legend,
                         lty = lty.legend,
                         pch = pch.legend,
                         col = col.legend,
                         pt.bg = pt.bg.legend,
                         bg="white",
                         cex=legend.cex[1],
                         pt.cex=legend.cex[2], 
                         ncol = legend.ncol)
    }
    if(add.arrow){
        n.decision <- length(ck)
        for(iDecision in 1:n.decision){
            graphics::arrows(x0=Idval[iDecision],y0=ylim[1],x1=Idval[iDecision],y1=yc[iDecision]-sep.arrow[1],lty=1,length=size.arrow,lwd=lwd.arrow)
            graphics::arrows(x0=Idval[iDecision],y0=yc[iDecision]+sep.arrow[2],x1=Idval[iDecision],y1=ylim[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            if(iDecision+1<=n.decision){
                graphics::arrows(x0=Idval[iDecision],y0=yc[iDecision]-sep.arrow[1],x1=Idval[iDecision+1],y1=ylim[1],lty=1,length=size.arrow,lwd=lwd.arrow)
                graphics::arrows(x0=Idval[iDecision+1],y0=ylim[2],x1=Idval[iDecision],y1=yc[iDecision]+sep.arrow[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            }else{
                graphics::arrows(x0=Idval[iDecision],y0=yc[iDecision]-sep.arrow[1],x1=Ival[iDecision+1],y1=ylim[1],lty=1,length=size.arrow,lwd=lwd.arrow)
                graphics::arrows(x0=Ival[iDecision+1],y0=ylim[2],x1=Idval[iDecision],y1=yc[iDecision]+sep.arrow[2],lty=1,length=size.arrow,lwd=lwd.arrow)
            }
        }
        graphics::arrows(x0=Ival[iDecision+1],y0=ylim[1],x1=Ival[iDecision+1],y1=ylim[2],lty=1,length=size.arrow,lwd=lwd.arrow)
    }
                                        # }}}    
}

