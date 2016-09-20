# Project: CSFA
# 
# Author: lucp8394
###############################################################################


#' Computing connectivity scores with Factor Analysis methodology.
#'
#' CSFA is a wrapper of multiple packages containing a factor analysis method.
#' These methods are used to derive the the connectivity scores of query gene signatures with one or multiple reference signatures.
#' CSFA will apply them, output the scores and immediately produce a number of meaningful plots interactively. The included methods are PCA and MFA from the \code{FactoMineR} package, FABIA from the \code{fabia} package and Sparse PCA/MFA from the \code{elasticnet} package.
#' Further, CSFA also contains an implementation of the Zhang and Gant score.
#' 
#' @references Abdi, H. et al. (2013), "Multiple factor analysis: principal component analysis for multitable and multiblock data sets," \emph{WIREs Comput Stat}, 1-31.
#' @references Hochreiter, S. et al., "FABIA: Factor Analysis for Bicluster acquisition," \emph{Bioinformatics}, 26, 1520-1527.
#' @references Lamb, J. et al. (2006), "The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease," \emph{Science}, 313, 1929-1934.
#' @references Zhang, S.-D. and Gant, T.W. (2008), "A simple and robust method for connecting small-molecule drugs using gene-expression signatures," \emph{BMC Bioinformatics}, 9, 10.
#' @references PAPER IN PROCESS: Connectivity Scores with Factor Analysis
#'
#' @docType package
#' @name CSFA
NULL