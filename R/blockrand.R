#' @title Block randomization
#' 
#' @param seed integer for the random seed generator state.
#' @param rand.block values within each block. The length of this vector determine the block size.
#' @param n number of observations.
#' 
#' @details Nothing
#' @return Vector or binary otcome.
#' @author Paul Blanche
#' 
#' @examples
#' blockrand(123,c(1,1,0,0),10)
#' @export
blockrand  <-  function(seed,rand.block,n){
    set.seed(seed)
    blocksize <- length(rand.block)
    block  <-  rep(1:ceiling(n/blocksize), each = blocksize)
    a1  <-  data.frame(block, rand=stats::runif(length(block)), envelope= 1: length(block))
    a2  <-  a1[order(a1$block,a1$rand),]
    a2$Z <- rep(rand.block,times = length(block)/blocksize)
    Z  <-  a2[order(a2$envelope),"Z"]
    return(Z[1:n])
}
