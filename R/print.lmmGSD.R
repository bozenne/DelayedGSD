## * print.lmmGSD
#' @export
print.lmmGSD <- function(x, ...){
    cat("       Analysis via the gls function \n \n ")
    cat("  Estimated treatment effect: ",x$name.coef," \n",sep = "")
    print(c("estimate" = x$delta$estimate, "se" = x$delta$se, "statistic" = x$delta$statistic, "df" = x$delta$df, "p.value" = x$delta$p.value))
    if(!is.null(x$information)){

        txt.legend <- "total     : all patients (including no observable outcome) \n"
        if(is.null(x$data.decision)){
            txt.legend <- c(txt.legend,
                            "decision  : patients with at least one observable outcome (including not yet observed values) \n")
            sample.size <- x$sample.size
        }else if(is.data.frame(x$data.decision)){
            txt.legend <- c("decision  : external dataset \n",
                            txt.legend)
            sample.size <- x$sample.size[c("decision","total","interim","interim.cc")]
        }else if(is.numeric(x$data.decision)){
            txt.legend <- c("decision  : user-input number of patients \n",
                            txt.legend)
            sample.size <- x$sample.size[c("decision","total","interim","interim.cc")]
        }
        
        cat("\n  Number of patients \n")
        print(sample.size)
        cat("\n  Estimated information \n")
        print(x$information)

        cat("\n",
            txt.legend,
            "interim   : patients with at least one observed outcome (full information) \n",
            "interim.cc: patients with at no missing outcome (complete case) \n")

    }
    return(NULL)
}
