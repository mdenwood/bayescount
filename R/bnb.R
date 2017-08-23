### Move to SuppDists

pbeta_nbinom <- function(q, a, c, d, ...){
  pghyper(q, -d, -a, c-1, ...)
}


pbnb1 <- function(q, k, alpha, beta, lower.tail = TRUE, log.p = FALSE){
  pghyper(q, -beta, -k, alpha-1, lower.tail = lower.tail, log.p = log.p)
}

pbnb2 <- function(q, k, alpha, beta, lower.tail = TRUE){
	
	stopifnot(length(q)==1)
	if(lower.tail){
		p <- .C(C_pbnb_lower_wrap, as.integer(q), as.double(k), as.double(alpha), as.double(beta), p = double(1))$p
	}
	return(p)
}

pghyper <-
function (q, a, k, N, lower.tail = TRUE, log.p = FALSE) 
{
    M <- max(length(q), length(a), length(k), length(N))
    q <- rep(q, length.out = M)
    a <- rep(a, length.out = M)
    k <- rep(k, length.out = M)
    N <- rep(N, length.out = M)
    if (lower.tail == TRUE) {
        value <- .C(C_pghyperR, as.integer(q), as.double(a), 
            as.double(k), as.double(N), as.integer(M), val = double(M))$val
    }
    else {
        value <- .C(C_ughyperR, as.integer(q), as.double(a), 
            as.double(k), as.double(N), as.integer(M), val = double(M))$val
    }
    if (log.p == TRUE) 
        value <- log(value)
    value
}

