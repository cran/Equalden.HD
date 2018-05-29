#' Rat data
#'
#' A microarray data set with 8038 logged gene expression levels measured on 5 rats. All rats were
#' subjected to the same treatment. The rows correspond to the genes and the columns refer to the rats.
#'
#' @usage data(Rat)
#'
#' @format A matrix with 8038 rows corresponding to the measured genes and 5 columns corresponding to the rats.
#'
#'
#' @references
#'
#' Davidson, L.A., Nguyen, D.V., Hokanson, R.M., Callaway, E.S., Isett, R.B., Turner, N.D., Dougherty, E.R., Wang, N., Lupton, J.R., Carroll, R.J., and Chapkin, R.S. (2004). Chemopreventive n-3
#' polyunsaturated fatty acids reprogram genetic signatures during colon cancer initiation and progression in the rat. Cancer Research, 64, 6797–6804.
#'
#' @examples
#' \dontshow{
#' data(Rat)
#' X <- Rat
#' k <- dim(X)[1]
#' s <- apply(X[1:100,], 1, density)
#' }
#' \donttest{
#' data(Rat)
#' X <- Rat
#' k <- dim(X)[1]
#' ### Estimated densities of logged gene expression levels.
#' s <- apply(X, 1, density)
#' ### Plot of estimated densities of 6 randomly selected genes.
#' set.seed(375)
#' rs <- sample(1:k, 6)
#' plot(s[[rs[1]]], main = "Kernel estimates for 6 randomly selected genes",
#'  xlab = "x", ylab = "density", xlim = c(-4, 2), ylim = c(0, 6))
#' for (i in 2:6){
#'   lines(s[[rs[i]]], col = i)
#' }
#' }

"Rat"
