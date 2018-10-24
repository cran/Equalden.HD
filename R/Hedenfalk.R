#' Hedenfalk data
#'
#' These data come from the breast cancer gene expression study of Hedenfalk et al. (2001). The data set
#' consists on 3226 gene expression levels measured on 7 patients with breast tumors having BRCA1
#' mutations (columns 1-7) and on 8 patients with breast tumors having BRCA2 mutations (columns 8-15).
#' The rows correspond to the genes and the columns refer to the patients.
#'
#' @usage data(Hedenfalk)
#'
#' @format A matrix with 3226 rows corresponding to the measured genes and 15 columns corresponding to the
#' patients. The first 7 columns contain the measures for the patients with BRCA1 mutations and the
#' patients with BRCA2 mutations are located in the last 8 columns.
#'
#'
#' @source \url{http://genomine.org/qvalue/results.txt}
#'
#' @references
#'
#' Hedenfalk, I., Duggan, D., Chen, Y., Radmacher, M, Bitter, M., Simon, R., Meltzer, P., Guterson, B., Esteller, M., Kallioniemi, O., Wilfond, B., Borg, A., Trent, J., Raffeld, M., Yakhini, Z., BenDor, A., Dougherty, E., Kononen, J., Buberdorf, L., Fehrle, W., Pittaluga, S., Gruvberger, G., Loman, N., Johannsson, O., Olsson, H., and Sauter, G.(2001), Gene-Expression Profiles in Hereditary Breast Cancer. New England Journal of Medicine 344 (8), 539-548.
#'
#' @examples
#' \dontshow{
#' data(Hedenfalk)
#' X <- Hedenfalk
#' k <- dim(X)[1]
#' s <- apply(X[1:100, 1:7], 1, density)
#' }
#'
#' \donttest{
#' data(Hedenfalk)
#' X <- Hedenfalk
#' k <- dim(X)[1]
#' ### Estimated densities of logged gene expression levels for patients with BRCA1 mutations.
#' s <- apply(X[, 1:7], 1, density)
#' ### Plot of estimated densities of 6 randomly selected genes.
#' set.seed(375)
#' rs <- sample(1:k, 6)
#' plot(s[[rs[1]]], main = "Kernel estimates for 6 randomly selected genes",
#'  xlab = "x", ylab = "density", xlim = c(-1, 4.3), ylim = c(0, 2))
#' for (i in 2:6){
#'   lines(s[[rs[i]]], col = i)
#' }
#' }
"Hedenfalk"
