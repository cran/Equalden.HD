h1 <- function(X, h) {
  p <- nrow(X)
  n <- ncol(X)
  Del <- X[, 1] - X[, 2:n]

  if (n > 2) {
    for (j in 2:(n - 1)) {
      Del <- cbind(Del, X[, j] - X[, (j + 1):n]) # we have to duplicated comparisons for this reason
      # in ans we multiply for 2.
    }
  }

  Del <- stats::dnorm(Del, sd = sqrt(2) * h)
  One <- matrix(1, n * (n - 1) / 2, 1) # n*(n-1)/2 número de comparacións sen contar as duplicadas.
  ans <- 2 * Del %*% One # we do the summation of the subscript j.
  ans <- ans / (n * (n - 1))
  as.vector(ans)
}


h3hat <- function(X, i, h) {
  p <- nrow(X)
  n <- ncol(X)
  Del <- rep(X[i, 1], len = p) - X # X_i1-x_kl para todo k,l.

  for (j in 2:n) {
    Del <- cbind(Del, rep(X[i, j], len = p) - X)
  }

  Del <- Del[(1:p)[(1:p) != i], ] # Sacamos as diferencias k=i.
  Del <- stats::dnorm(Del, sd = sqrt(2) * h)
  sum(Del) / (n ^ 2 * (p - 1))
}



teststat <- function(h, X) {
  p <- nrow(X)
  n <- ncol(X)

  h1vec <- 1:p
  h3est <- rep(0, len = p)
  h1vec <- h1(X, h)

  for (j in 1:p) {
    h3est[j] <- h3hat(X, j, h)
  }

  SW <- mean(h1vec)
  SB <- mean(h3est)
  list(SW - SB)
}


hseu <- function(X, h, es) {
  p <- nrow(X)
  n <- ncol(X)

  h1vec <- 1:p
  h3est <- rep(0, len = p)
  h1vec <- h1(X, h)

  for (j in 1:p) {
    h3est[j] <- h3hat(X, j, h)
  }

  sum1 <- rep(0, p)

  for (i in 1:p) {
    sum1[i] <- sum(h1vec[-i]) * (1 / (p - 1)) * 0.5
  }

  sum <- sum1 + 0.5 * h1vec - h3est
  return(sum - es)
}



bOptU <- function(influ, weights = c("parzen", "bartlett")) {
  weights <- match.arg(weights)
  n <- length(influ)

  ## Parameters for adapting the approach of Politis and White (2004)
  kn <- max(5, ceiling(log10(n)))
  lagmax <- ceiling(sqrt(n)) + kn

  ## Determine L
  L <- Lval(matrix(influ), method = min)

  ## Compute gamma.n
  gamma.n <- as.numeric(stats::ccf(influ, influ, lag.max = L,
                                   type = "covariance", plot = FALSE)$acf)

  sqrderiv <- switch(weights,
                     bartlett = 143.9977845,
                     parzen = 495.136227)
  integralsqrker <- switch(weights,
                           bartlett = 0.5392857143,
                           parzen = 0.3723388234)

  ft <- flattop(-L:L / L)
  Gamma.n.2 <- sqrderiv / 4 * sum(ft * (-L:L) ^ 2 * gamma.n) ^ 2
  Delta.n <- integralsqrker * 2 * sum(ft * gamma.n) ^ 2
  ln.opt <- (4 * Gamma.n.2 / Delta.n * n) ^ (1 / 5)
  round(max(ln.opt, 1))
}



### The next functions are used to compute the function \varphi in the
### variance estimator, equation (20), which is defined in the first equation in page 4.

parzen <- function(x) {
  ifelse(abs(x) <= 1/2, 1 - 6 * x^2 + 6 * abs(x)^3,
         ifelse(1/2 <= abs(x) & abs(x) <= 1, 2 * (1 - abs(x))^3, 0))
}

pdfsumunif <- function(x,n) {
  nx <- length(x)

  .C("pdf_sum_unif",
     as.integer(n),
     as.double(x),
     as.integer(nx),
     pdf = double(nx),
     PACKAGE = "Equalden.HD")$pdf
}


convrect <- function(x, n) {
  pdfsumunif(x + n/2, n) / pdfsumunif(n / 2, n)
}


flattop <- function(x, a = 0.5) {
  pmin(pmax((1 - abs(x)) / (1 - a), 0), 1)
}


sigmaf <- function(hseudo, ln_opt) {
  phi <- {
    function(x) convrect(x * 4, 8)
  }
  sum <- 0
  p <- length(hseudo)

  for (i in 1:p) {

    for (j in 1:p) {

      sum <- sum + phi((i - j) / ln_opt) * hseudo[i] * hseudo[j]
    }
  }

  sum / p
}



## Adapted from Matlab code by A. Patton and the R translation
## by C. Parmeter and J. Racine
mval <- function(rho, lagmax, kn, rho.crit) {
  ## Compute the number of insignificant runs following each rho(k),
  ## k=1,...,lagmax.
  num.ins <- sapply(1:(lagmax-kn+1),
                    function(j) sum((abs(rho) < rho.crit)[j:(j+kn-1)]))

  ## If there are any values of rho(k) for which the kn proceeding
  ## values of rho(k+j), j=1,...,kn are all insignificant, take the
  ## smallest rho(k) such that this holds (see footnote c of
  ## Politis and White for further details).
  if(any(num.ins == kn)) {
    return(which(num.ins == kn)[1])
  } else {
    ## If no runs of length kn are insignificant, take the smallest
    ## value of rho(k) that is significant.
    if(any(abs(rho) > rho.crit)) {
      lag.sig <- which(abs(rho) > rho.crit)
      k.sig <- length(lag.sig)

      if(k.sig == 1)
        ## When only one lag is significant, mhat is the sole
        ## significant rho(k).
        return(lag.sig)
      else
        ## If there are more than one significant lags but no runs
        ## of length kn, take the largest value of rho(k) that is
        ## significant.
        return(max(lag.sig))
    }
    else
      ## When there are no significant lags, mhat must be the
      ## smallest positive integer (footnote c), hence mhat is set
      ## to one.
      return( 1 )
  }
}

Lval <- function(x, method = mean) {
  x <- matrix(x)
  n <- nrow(x)
  d <- ncol(x)

  ## parameters for adapting the approach of Politis and White (2004)
  kn <- max(5, ceiling(log10(n)))
  lagmax <- ceiling(sqrt(n)) + kn
  rho.crit <- 1.96 * sqrt(log10(n) / n)

  m <- numeric(d)

  for (i in 1:d) {
    rho <- stats::acf(x[, i], lag.max = lagmax, type = "correlation",
                      plot = FALSE)$acf[-1]
    m[i] <- mval(rho, lagmax, kn, rho.crit)
  }
  return(2 * method(m))
}

covariance <- function(hseudo, m) {
  p <- length(hseudo)
  c <- rep(0, m)

  for (i in 1:m) {
    sum <- 0

    for (j in 1:(p-i)){
      sum <- sum + hseudo[j] * hseudo[j + i]
    }

    c[i] <- sum

  }
  return((1/p)*c)
}

variance <- function(hseudo) {
  p <- length(hseudo)
  var <- 0
  for (j in 1:(p)){
    var <- var + hseudo[j] * hseudo[j]
  }

  return((1 / p) * var)
}

statistic <- function(c, hseudo) {
  m <- length(c)
  part2 <- 0

  for (i in 1:m){ #dyn.load
    part2 <- part2 + (1- (i/(m+1)))*c[i]
  }

  statistic <- variance(hseudo) + 2 * part2
  return(statistic)
}

teststat2 <- function(h, X) {
  p <- nrow(X)
  n <- ncol(X)
  h1vec <- 1:p
  h3est <- rep(0, len = p)
  h1vec <- h1(X, h)

  for (j in 1:p) {
    h3est[j] <- h3hat(X, j, h)
  }

  SW <- mean(h1vec)
  SB <- mean(h3est)
  sig2 <- stats::var(h1vec - 2 * h3est)
  list(sqrt(p) * (SW - SB) / sqrt(sig2), SW - SB, sig2 = sig2) # quitar sig2
  # return(list(sig2 = sig2))
}



h1 <- function(X, h) {
  p <- nrow(X)
  n <- ncol(X)
  Del <- X[, 1] - X[, 2:n]
  if (n > 2) {
    for (j in 2:(n - 1)) {
      Del <- cbind(Del, X[, j] - X[, (j + 1):n]) # we have to duplicated comparisons for this reason
      # in ans we multiply for 2.
    }
  }
  Del <- stats::dnorm(Del, sd = sqrt(2) * h)
  One <- matrix(1, n * (n - 1) / 2, 1) # n*(n-1)/2 n?mero de comparaci?ns sen contar as duplicadas.
  ans <- 2 * Del %*% One # we do the summation of the subscript j.
  ans <- ans / (n * (n - 1))
  as.vector(ans)
}

### Function to compute the function h2 in (1), more precisely, for each this
### function reports (1/(p-1))*\sum_{k=1,k\not=i} h_2(X_i,Xk).


h3hat <- function(X, i, h) {
  p <- nrow(X)
  n <- ncol(X)
  Del <- rep(X[i, 1], len = p) - X # X_i1-x_kl para todo k,l.
  for (j in 2:n) {
    Del <- cbind(Del, rep(X[i, j], len = p) - X)
  }
  Del <- Del[(1:p)[(1:p) != i], ] # sacamos as diferencias k=i.
  Del <- stats::dnorm(Del, sd = sqrt(2) * h)
  sum(Del) / (n ^ 2 * (p - 1))
}


### Using the previous function, the next function computes the statistic S


teststat <- function(h, X) {
  p <- nrow(X)
  n <- ncol(X)
  h1vec <- 1:p
  h3est <- rep(0, len = p)
  h1vec <- h1(X, h)
  for (j in 1:p) {
    h3est[j] <- h3hat(X, j, h)
  }
  SW <- (h1vec)
  SB <- (h3est)
  list(SW - SB)
}


