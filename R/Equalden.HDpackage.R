#' Package ‘Equalden.HD’
#'
#' Documentation for package ‘Equalden.HD’ version 1.0
#'
#' @description
#' This package implements three different methods to test the null hypothesis that a large number k of
#' samples have a common density. The sample size can be as small as 2. These methods are particularly
#' well suited to the low sample size, high dimensional setting (n << k). The first method, proposed by
#' Zhan and Hart (2012), was developed to test the null hypothesis when the samples are independent of
#' each other. The other tests, proposed by Cousido-Rocha et al. (2018), are adaptations of the test in Zhan
#' and Hart (2012) for the setting in which the samples are weakly dependent. The standarized version of
#' each test statistic and its p-value are computed among other things.
#'
#'
#' @details
#' \itemize{
#' \item{Package: Equalden.HD}
#' \item{Version: 1.1}
#' \item{Maintainer: Marta Cousido Rocha \email{martacousido@@uvigo.es}}
#' \item{License: GPL-2}
#' }
#'
#' @return
#' \itemize{
#' \item{Equalden.test.HD: Performs the k-sample test proposed in Zhan and Hart (2012) for the setting of low sample size, large
#' dimension and independent samples, and its adaptions to dependent samples proposed in Cousido-Rocha
#' et. al (2018).}
#' }
#'
#' @author
#' \itemize{
#' \item{Cousido Rocha, Marta.}
#' \item{Soage González, José Carlos.}
#' \item{de Uña-Álvarez, Jacobo.}
#' \item{D. Hart, Jeffrey.}
#' }
#'
#' @section Acknowledgements:
#' This work has received financial support of the Call 2015 Grants for
#' PhD contracts for training of doctors of the Ministry of Economy and Competitiveness,
#' co-financed by the European Social Fund (Ref. BES-2015-074958).
#' The authors acknowledge support from  MTM2014-55966-P project, Ministry of Economy and
#' Competitiveness, and MTM2017-89422-P project, Ministry of Economy, Industry
#' and Competitiveness, State Research Agency, and Regional Development Fund, UE.
#' The authors also acknowledge the financial support provided by the SiDOR research group
#' through the grant Competitive Reference Group, 2016-2019 (ED431C 2016/040),
#' funded by the ``Consellería de Cultura, Educación e Ordenación Universitaria.
#' Xunta de Galicia''. José Carlos Soage was supported by Red Tecnológica de
#' Matemática Industrial (Red TMATI), Cons. de Cultura, Educación e OU, Xunta de Galicia
#' (ED341D R2016/051) and by Grupos de Referencia Competitiva, Consolidación y
#' Estructuración de Unidades de Investigación Competitivas del SUG, Cons. de Cultura,
#' Educación e OU, Xunta de Galicia (GRC ED431C 2016/040).
#'
#' @references
#' \itemize{
#' \item{Cousido-Rocha, M., de Uña-Álvarez, J., and Hart, J.(2018). Testing equality of a large number of densities under mixing conditions. Preprint.}
#' \item{Zhan, D., Hart, J. (2012). Testing equality of a large number of densities. Biometrika, 99, 1-17.}
#' }
#'
#'
"_PACKAGE"
#> [1] "_PACKAGE"
