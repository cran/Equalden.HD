#' @title A test for the equality of a high dimensional set of densities
#'
#' @description
#' Performs the k-sample test proposed by Zhan and Hart (2012) for the low sample size, high dimensional
#' setting with independent samples, and its extensions for dependent samples proposed by Cousido-Rocha
#' et al. (2018).
#'
#' @usage
#' Equalden.test.HD(X, method = c("indep", "dep.boot", "dep.spect"))
#'
#' @param X A matrix where each row is one of the k-samples.
#' @param method the k-sample test. By default the "dep.spect" method is computed. See details.
#'
#' @details
#' The function implements the k-sample test proposed by Zhan and Hart (2012), method="indep", and
#' its extensions for dependent data proposed by Cousido-Rocha et al. (2018), method="dep.boot" and
#' "dep.spect". The method proposed by Zhan and Hart (2012) serves to test the null hypothesis that the
#' k-samples have a common distribution. It is suitable when the k samples are independent and the number
#' of samples k is large, and it works for sample sizes as small as 2. The statistic in Zhan and Hart (2012)
#' is based on a comparison between the k sample-specific kernel density estimates and the kernel density
#' estimate computed from the pooled sample. An alternative expression of this statistic shows that it can
#' be interpreted as a difference between the intra-samples variability and the inter-samples variability.
#' This statistic is standarized using a variance estimator which is valid for independent samples. The asymptotic
#' normality (when k tends to infinity) of the standardized version of the statistic is used to compute the
#' corresponding p-value. Cousido-Rocha et al. (2018) proposed two corrections of the test of Zhan and
#' Hart (2012) for dependent samples. These tests standarize the statistic proposed in Zhan and Hart (2012)
#' by using variance estimators which are suitable when the samples are weakly dependent. The method
#' "dep.boot" implements the dependent multiplier bootstrap to estimate the variance, whereas the method
#' "dep.spect" uses a variance estimator based on the spectral analysis theory. Both tests perform similarly,
#' but the "dep.spect" test tends to be computationally more efficient than the "dep.boot" test. Cousido-
#' Rocha et al. (2018) showed through simulations that, for independent samples, the tests "dep.boot" and
#' "dep.spect" may be more powerful than the test in Zhan and Hart (2012) despite of being protected
#' against possible dependences.
#' On the other hand, the statistic can be written as a sum of k individual statistics, each of them
#' measures the difference between the intra-sample variability of the corresponding variable and the inter-
#' samples variability. Whether the null hypothesis is rejected, an exploratory analysis of such individual
#' statistics can help to guess which genes are not equally distributed.
#'
#' @return A list containing the following components:
#' \item{standarized statistic: }{the value of the standarized statistic.}
#' \item{p.value: }{the p-value for the test.}
#' \item{statistic: }{the value of the statistic.}
#' \item{variance: }{the value of the variance estimator.}
#' \item{m: }{number of significant lags for the variance estimator if the method is "dep.spect" or "dep.boot". Null if the method is "indep" since no correction for dependence is required in this case.}
#' \item{k: }{number of samples or populations.}
#' \item{n: }{sample size.}
#' \item{method: }{a character string indicating what k-test was performed.}
#' \item{I.statistics: }{the k individual statistics.}
#' \item{data.name: }{a character string giving the name of the data.}
#'
#'
#' @author
#' \itemize{
#' \item{Marta Cousido-Rocha}
#' \item{José Carlos Soage González}
#' \item{Jacobo de Uña-Álvarez}
#' \item{D. Hart, Jeffrey.}
#' }
#'
#' @references
#' \itemize{
#' \item{Cousido-Rocha, M., de Uña-Álvarez, J., and Hart, J.(2018). Testing equality of a large number of densities under mixing conditions. Preprint.}
#' \item{Zhan, D., Hart, J. (2012). Testing equality of a large number of densities. Biometrika, 99, 1-17.}
#' }
#'
#' @examples
#' \dontshow{
#' n <- 2
#' k <- 100
#' set.seed(1234)
#' X <- matrix(rnorm(n * k), ncol = 2)
#' res <- Equalden.test.HD(X,  method = "indep")
#' res
#'
#' ### The statistic and the variance estimator
#' res$statistic
#' res$variance
#' }
#' \donttest{
#' n <- 2
#' k <- 100
#' set.seed(1234)
#' X <- matrix(rnorm(n * k), ncol = 2)
#' res <- Equalden.test.HD(X,  method = "indep")
#'
#' res
#' ### The statistic and the variance estimator
#' res$statistic
#' res$variance
#' ### The number of samples and sample size
#' res$k
#' res$n
#'
#' ### Real data analysis. We test the null hypothesis that 1000 randomly selected genes
#' ### measured in patients with BRCA2 mutations have a common distribution. We use the test
#' ### proposed in Cousido-Rocha et al. (2018) since correlation among expression levels of
#' ### different genes on the same individual is expected.
#' data(Hedenfalk)
#' X <- Hedenfalk
#' k <- dim(X)[1]
#' ### We eliminate the additive patients effects by substracting to each column its sample mean.
#' BRCA2 <- sweep(X[, 8:15], 2, apply(X[, 8:15], 2, mean))
#' set.seed (1234)
#' se<-1000
#' ind <- sample(1:k, se)
#' res1 <- Equalden.test.HD(BRCA2[ind, ], method = "dep.boot")
#' res1
#' res2 <- Equalden.test.HD(BRCA2[ind, ], method = "dep.spect")
#' res2
#' ### The null hypothesis is rejected using both methods. Then we plot the individual statistics
#' ### and highlight the 100 most extreme values.
#' cu <- 100
#' I.statistics.sorted <- sort(res1$I.statistics)
#' cv <- I.statistics.sorted[se-cu+1]
#' ind2 <- which(res1$I.statistics >= cv)
#' plot(1:se, res1$I.statistics, xlim = c(0, se), ylim = c(min(res1$I.statistics),
#'                                                         max(res1$I.statistics)),
#'      xlab = "Genes", ylab = "statistic", main = "Individual statistics")
#' points(ind2, res1$I.statistics[ind2], col = "red")
#' ### We zoom the plot in the following way since some individual statistics report extreme
#' ### negative values in this data.
#' plot(1:se, res1$I.statistics, xlim = c(0, se), ylim = c(0, max(res1$I.statistics)),
#'      xlab = "Genes", ylab = "statistic", main = "Individual statistics")
#' points(ind2, res1$I.statistics[ind2], col = "red")
#' }
#' @importFrom  stats var
#' @useDynLib Equalden.HD, .registration = TRUE
#' @export
Equalden.test.HD <- function(X, method = c("indep", "dep.boot", "dep.spect")) {
  cat("Call:", "\n")
  print(match.call())
  if(missing(method)) {
    method <- "dep.spect"
    cat("'dep.spect' method used by default")
  }
  method <- match.arg(method)
  DNAME <- deparse(substitute(X))
  METHOD <- "A test for the equality of a high dimensional set of densities"

  match.arg(method)

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


  n <- ncol(X)

  if(n < 2){
    stop("Number of columns must be at least 2.")
  }

  p <- nrow(X)

  if(p < 2){
    stop("Number of rows must be at least 2.")
  }


  #=============================================================================
  # Sp
  #=============================================================================

  ### The next lines compute the bandwidht, see equation (23).
  c1 <- 1 / p
  c2 <- 1.144 * (n ^ (-1 / 5))
  si <- apply(X, 1, stats::var)
  spool <- sqrt(c1 * sum(si))
  h <- spool * c2

  ### The next line computes the statistic S
  e <- unlist(teststat(h, X))

  ### The next lines computes the pseudo-observations, the parameter
  ### \widehat{l}_p^{opt} and the variance estimator, equation (20).
  hseudo <- hseu(X, h, e)

  ln_opt <- bOptU(hseudo)
  sigma <- sigmaf(hseudo, ln_opt)

  ### Then we compute S_p, the standarized version of our test statistic
  eso <- (sqrt(p) * (e)) / (2 * sqrt((sigma)))

  ### Then we compute the p-value using the asymptotic normality.
  pvalor1 <- 1 - stats::pnorm(eso)


  #=============================================================================
  # \tilde{S}_p
  #=============================================================================

  ### The next line computes the value of the parameter m, which is assumed
  ### equal to L_p in the expression of \widehat{l}_p^{opt}
  m <- Lval(hseudo, method = min)

  ### The next lines computes the estimator of the covariance function in (22)
  c <- covariance(hseudo, m)
  ### Then, we compute the variance estimator in (22).
  sigma2 <- statistic(c, hseudo)

  ### Then we compute S_p, the standarized version of our test statistic
  esa <- (sqrt(p) * (e)) / (2 * sqrt((sigma2)))

  ### Then we compute the p-value using the asymptotic normality.
  pvalor2 <- 1 - stats::pnorm(esa)


  #=============================================================================
  # ZH test
  #=============================================================================

  ### The next line computes the standarized version of ZH test.
  a <- teststat2(h, X)
  s <- a[1]
  sig2 <- a$sig2

  ### Then we compute the p-value using the asymptotic normality.
  pvalor3 <- 1 - stats::pnorm(unlist(s))


  #=============================================================================
  # I.statistics
  #=============================================================================

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

  c1 <- 1 / p
  c2 <- 1.144 * n ^ (-1 / 5)
  si <- apply(X, 1, var)
  spool <- sqrt(c1 * sum(si))
  h <- spool * c2


  stats_ind <- teststat(h, X) # I.statistics is equal to stats_ind



  statistic <- switch(method, dep.boot = eso, dep.spect = esa, indep = unlist(s))
  names(statistic) <- "standarized statistic"

  statistic2 <- switch(method, dep.boot = eso, dep.spect = esa, indep = unlist(s))

  p.value <- switch(method, dep.boot = pvalor1, dep.spect = pvalor2, indep = pvalor3)

  met <- switch(method, dep.boot = "dep.boot", dep.spect = "dep.spect", indep = "indep")

  variance <- switch(method, dep.boot = sigma, dep.spect = sigma2, indep = sig2)

  m <- switch(method, dep.boot = ln_opt, dep.spect = m, indep = NULL)

  RVAL <- list(statistic = statistic, p.value = p.value, method = METHOD,
               data.name = DNAME, sample.size = n, method1 = met)

  RVAL2 <- list(standarized.statistic = statistic2, p.value = p.value,
                statistic = e, variance = variance, m = m, k = p,
                n = n, method = met, data.name = DNAME, I.statistics = unlist(stats_ind))
  class(RVAL) <- "htest"

  print(RVAL)
  return(invisible(RVAL2))
}

