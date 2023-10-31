#' @title rho-family spending functions (Kim-DeMets) for alpha and beta
#'
#' @param I Information at analysis to be evaluated
#' @param rho parameter for the spending function
#' @param beta_or_alpha total alpha or beta to be spent
#' @param Info.max Maximum information needed at end of trial
#'
#' @export
ErrorSpend <- function(I,
                       rho,
                       beta_or_alpha,
                       Info.max){
    if(length(I)==0){
        stop("Information at analysis is missing.")
    }
    return(beta_or_alpha*min(1,(I/Info.max)^rho))
}
