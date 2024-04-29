### plot.gridDelayedGSD.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 23 2024 (09:09) 
## Version: 
## Last-Updated: apr 24 2024 (09:08) 
##           By: Brice Ozenne
##     Update #: 15
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * plot.gridDelayedGSD
##' @title Display p-values
##' @description Display p-value over the sample space.
##'
##' @param object outpout of gridDelayedGSD.
##' @param col [character vector] color used to display the arrows.
##' @param xlab [character] label for the x-axis.
##' @param ylab [character] label for the y-axis.
##' @param title [character] title of the plot.
##' @param xlim [numeric vector of length 2] minimum and maximum value on the x-axis.
##' @param ylim [numeric vector of length 2] minimum and maximum value on the y-axis.
##' @param pch [integer, >0] shape used to display the boundaries at decision and final.
##' @param cex [numeric, >0] size of the shape used to display the boundaries at decision and final.
##' @param cex.lab [numeric, >0] size of the x-axis and y-axis labels.
##' @param cex.p [numeric, >0] size of the used to display the p-values.
##' @param font.p [numeric, >0] font used to display the p-values.
##' @param digits [integer, >0] number of digit to be used when displaying the p-values.
##' 
##' 
##' @export
plot.gridDelayedGSD <- function(object,
                                col = NULL,
                                xlab = "Stage",
                                ylab = "Test statistic",
                                title = NULL,
                                xlim = NULL,
                                pch = 20,
                                cex = 2,
                                cex.lab = 1.2,
                                cex.p = 1,
                                font.p = 2,
                                lwd = 2,
                                digits = 2){

    ## ** graphical display
    kMax <- max(object$k)
    if(is.null(xlim)){
        xlim <- c(1-0.1*digits,kMax+0.1*digits)
    }
    n.col <- 1 + 4*(kMax-1)
    if(is.null(col)){
        col <- hcl.colors("Blues", n = n.col+3)
    }else if(length(col)!=n.col){
        stop("Argument color should have length ",n.col,". \n")
    }
    index.change0 <- which(diff(object$k)!=0)
    index.boundary <- c(index.change0[1:(kMax-1)],which.min(abs(object$p.value-0.025)))
    object$p.display <- paste0(round(100*object$p.value, digits = digits),"%")


    plot(object$k,object$z, type = "l", xlab = xlab, ylab = ylab, cex.lab=cex.lab,
         axes = FALSE, xlim = xlim, main = title)
    for(iStage in 1:(kMax-1)){
        if(object[index.boundary[iStage],"fixC"]){
            points(object$k[index.boundary[iStage]-1],object$z[index.boundary[iStage]-1], pch = pch, cex = cex)
            points(object$k[index.boundary[iStage]],object$z[index.boundary[iStage]], pch = pch, cex = cex)
        }else{
            points(object$k[index.boundary[iStage]],object$z[index.boundary[iStage]], pch = pch, cex = cex)
        }
    }
    points(object$k[index.boundary[kMax]],object$z[index.boundary[kMax]], pch = pch, cex = cex)
    axis(1, at = 0:kMax)
    axis(2)
    for(iStage in 1:(kMax-1)){
        if(object[index.boundary[iStage],"fixC"]){
            arrows(x0 = object[c(1,index.change0+1)[iStage],"k"], y0 = object[c(1,index.change0+1)[iStage],"z"], 
                   x1 = object[index.change0[iStage]-1,"k"], y1 = object[index.change0[iStage]-1,"z"],
                   col = col[1+2*(iStage-1)], lwd = lwd)
            points(x = object[c(index.change0[iStage]-1, index.change0[iStage]),"k"], y = object[c(index.change0[iStage]-1, index.change0[iStage]),"z"],
                   col = "red", lwd = lwd, type = "l")
        }else{
            arrows(x0 = object[c(1,index.change0+1)[iStage],"k"], y0 = object[c(1,index.change0+1)[iStage],"z"], 
                   x1 = object[index.change0[iStage],"k"], y1 = object[index.change0[iStage],"z"],
                   col = col[1+2*(iStage-1)], lwd = lwd)
        }
        arrows(x0 = object[index.change0[iStage],"k"], y0 = object[index.change0[iStage],"z"],
               x1 = object[index.change0[iStage]+1,"k"], y1 = object[index.change0[iStage]+1,"z"],
               col = col[2+2*(iStage-1)], lwd = lwd)
        iSeq <- seq(c(1,index.change0+1)[iStage],index.change0[iStage]-1,by=1)
        text(x = object[iSeq,"k"], y = object[iSeq,"z"],label = object[iSeq,"p.display"], pos = 2, cex = cex.p, font = font.p)
        text(x = object[index.change0[iStage],"k"], y = 0.35*object[max(iSeq),"z"]+0.65*object[index.change0[iStage],"z"],
             label = object[index.change0[iStage],"p.display"], pos = 2, col = ifelse(object[index.boundary[iStage],"fixC"],"red","black"), cex = cex.p, font = font.p)
    }
    arrows(x0 = object[index.change0[kMax-1]+1,"k"], y0 = object[index.change0[kMax-1]+1,"z"],
           x1 = object[index.change0[kMax],"k"], y1 = object[index.change0[kMax],"z"],
           col = col[1+2*(kMax-1)], lwd = lwd)
    text(x = object[(index.boundary[kMax-1]+1):index.change0[kMax],"k"], y = object[(index.boundary[kMax-1]+1):index.change0[kMax],"z"],
         label = object[(index.boundary[kMax-1]+1):index.change0[kMax],"p.display"], pos = 4, cex = cex.p, font = font.p)

    for(iStage in 1:(kMax-1)){ ## iStage <- 1
        arrows(x0 = object[index.change0[kMax+iStage-1],"k"], y0 = object[index.change0[kMax+iStage-1],"z"],
               x1 = object[index.change0[kMax+iStage-1]+1,"k"], y1 = object[index.change0[kMax+iStage-1]+1,"z"],
               col = col[1+2*(kMax+iStage-1)], lwd = lwd)
        arrows(x0 = object[index.change0[kMax+iStage-1]+1,"k"], y0 = object[index.change0[kMax+iStage-1]+1,"z"],
               x1 = object[c(index.change0,NROW(object))[kMax+iStage],"k"], y1 = object[c(index.change0,NROW(object))[kMax+iStage],"z"],
               col = col[2+2*(kMax+iStage-1)], lwd = lwd)
        text(x = object[index.change0[kMax+iStage-1]+1,"k"], y = object[index.change0[kMax+iStage-1]+1,"z"]*1.05,
             label = object[index.change0[kMax+iStage-1]+1,"p.display"], pos = 2, cex = cex.p, font = font.p)
        iSeq <- seq(index.change0[kMax+iStage-1]+2,c(index.change0,NROW(object))[kMax+iStage],by=1)
        text(x = object[iSeq,"k"], y = object[iSeq,"z"], label = object[iSeq,"p.display"], pos = 2, cex = cex.p, font = font.p)
    }

    return(invisible(NULL))
}



##----------------------------------------------------------------------
### plot.gridDelayedGSD.R ends here
