
#' Duplicate gene pairs in the genome of Arabidopsis thaliana
#'
#' Duplicate gene pairs were identified and classified by duplication
#' mode with the \strong{doubletrouble} package.
#'
#' @name ath_dups
#' @format A data frame with the following variables:
#' \itemize{
#'   \item{\emph{dup1} Character, duplicated gene 1.}
#'   \item{\emph{dup2} Character, duplicated gene 2.}
#'   \item{\emph{type} Factor, duplication mode.}
#' }
#' @examples
#' data(ath_dups)
#' @usage data(ath_dups)
"ath_dups"


#' Expression data for experiment SRP213876
#'
#' Expression data (in TPM) and sample metadata were obtained from
#' Refine Bio.
#'
#' @name ath_se
#' @format A `SummarizedExperiment` object with an assay named \strong{TPM}
#' and three colData columns:
#' \itemize{
#'   \item{\emph{experiment} Character, experiment accession.}
#'   \item{\emph{tissue} Character, tissue name.}
#'   \item{\emph{treatment} Character, treatment.}
#' }
#' @examples
#' data(ath_se)
#' @usage data(ath_se)
"ath_se"


#' Genes and modules for a coexpression network inferred from SRP201971
#'
#' Expression data (in TPM) and sample metadata were obtained from
#' Refine Bio. The signed gene coexpression network was inferred with
#' the \strong{BioNERO} package.
#'
#' @name genes_and_modules
#' @format A 2-column data frame with the following variables:
#' \itemize{
#'   \item{\emph{Genes} Character, gene ID.}
#'   \item{\emph{Modules} Character, module name.}
#' }
#' @examples
#' data(genes_and_modules)
#' @usage data(genes_and_modules)
"genes_and_modules"


#' Module eigengenes for a coexpression network inferred from SRP201971
#'
#' Expression data (in TPM) and sample metadata were obtained from
#' Refine Bio. The signed gene coexpression network was inferred with
#' the \strong{BioNERO} package.
#'
#' @name ME
#' @format A data frame with samples in rows and module eigengenes in columns.
#' @examples
#' data(ME)
#' @usage data(ME)
"ME"

