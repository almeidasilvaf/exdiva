

#' Calculate correlation coefficients and permutation P-values for gene pairs
#'
#' @param se A `SummarizedExperiment` object with expression data.
#' @param dups A data frame of duplicate pairs, with duplicated genes
#' for each pair in columns 1 and 2. Additional columns are allowed, and will
#' be ignored if present.
#' @param method Character indicating which correlation coefficient to compute.
#' One of 'pearson', 'spearman', or 'kendall'. Default: 'spearman'.
#' 
#' @return A data frame exactly as the one passed in \strong{dups},
#' but with an extra column named `cor` (numeric) with correlation 
#' coefficients. Pairs will have NA in this column if i. at least of the
#' genes in the pair is not present in the assay slot of \strong{se}; or ii.
#' at least of the genes in the pair has a standard deviation of 0.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor sd
#' @rdname calc_cor
#' @examples
#' # Load data
#' data(ath_se)
#' data(ath_dups)
#' 
#' # Subset data
#' dups <- ath_dups[ath_dups$type == "TD", ]
#' se <- ath_se
#' 
#' # Calculate correlation
#' cor_df <- calculate_cor(se, dups, method = "pearson")
calculate_cor <- function(se, dups, method = "spearman") {
    
    exp <- assay(se)
    genes <- rownames(se)
    
    # Remove genes with zero standard deviation
    zero_sd <- apply(exp, 1, sd)
    zero_sd <- names(zero_sd[zero_sd == 0])
    if(length(zero_sd) > 0) { 
        genes <- genes[!genes %in% zero_sd]
    }
    
    # Calculate correlation between gene pairs
    cor_df <- Reduce(rbind, lapply(seq_len(nrow(dups)), function(x) {
        gene1 <- dups[x, 1]
        gene2 <- dups[x, 2]
        
        corr <- NA
        if(all(c(gene1, gene2) %in% genes)) {
            corr <- cor(exp[gene1, ], exp[gene2, ], method = method)
        }
        
        return(cbind(dups[x, ], cor = corr))
    }))
    
    return(cor_df)
}


#' Calculate \eqn{\tau} (Tau) index of tissue specificity
#' 
#' @param aggregated_exp A gene expression matrix with genes in rows
#' and tissues in columns as returned by \code{aggregate_to_tissue()}.
#' @param log Logical indicating whether the matrix of aggregated
#' gene expression per tissue should be log-transformed
#' using \eqn{log_2(x+1)} before calculating \eqn{\tau}.
#' 
#' @return A data frame with the following variables:
#' \itemize{
#'   \item{\emph{gene} Character, gene ID.}
#'   \item{\emph{tau} Numeric, \eqn{\tau} index of tissue specificity.}
#'   \item{\emph{class} Factor, expression-based class. One of 'Null', 
#'   'Weak', 'Broad', 'Specific' (see Details field).}
#' }.
#' 
#' @details
#' Besides calculating \eqn{\tau} indices, this function classifies genes into
#' expression-based classes, namely:
#' \itemize{
#'   \item{\strong{Null} Expression <1 in all tissues.}
#'   \item{\strong{Weak} Expression <5 in all tissues.}
#'   \item{\strong{Broad} \eqn{\tau} <0.85.}
#'   \item{\strong{Specific} \eqn{\tau} >=0.85.}
#' }
#'
#' @rdname calculate_tau
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats median
#' @examples
#' # Load data
#' data(ath_se)
#' se <- ath_se
#' 
#' # Simulate different tissues
#' SummarizedExperiment::colData(se)$tissue <- c(rep("seed", 3), rep("leaf", 3))
#'
#' # Aggregate expression to tissue level
#' aexp <- aggregate_to_tissue(se, "tissue")
#' 
#' # Calculate Tau
#' tau_df <- calculate_tau(aexp)
calculate_tau <- function(aggregated_exp, log = TRUE) {
    
    fexp <- aggregated_exp
    if(log) { fexp <- log2(aggregated_exp + 1)}
    
    # Calculate Tau
    tau <- apply(fexp, 1, function(x) {
        
        v <- NA
        if(all(!is.na(x)) & min(x, na.rm = TRUE) >= 0) {
            v <- 0
            if(max(x) != 0) {
                x <- (1-(x / max(x, na.rm = TRUE)))
                v <- sum(x, na.rm = TRUE)
                v <- v / (length(x) - 1)
            }
        }
        
        return(v)
    })
    tau_df <- data.frame(gene = names(tau), tau = tau, row.names = NULL)
    
    # Classify genes into expression-based categories
    tau_df <- classify_tau(aggregated_exp, tau_df)

    return(tau_df)
}


#' Compare \eqn{\tau} indices for gene pairs
#'
#' @param tau_df A 3-column data frame with \eqn{\tau} indices and expression
#' classes for each gene as returned by \code{calculate_tau()}.
#' @param dups A data frame of duplicate pairs, with duplicated genes
#' for each pair in columns 1 and 2. Additional columns are allowed, and will
#' be ignored if present.
#'
#' @return A data frame with the same columns as the input to \strong{dups},
#' but with the following extra columns:
#' \itemize{
#'   \item{\emph{tau_dup1} Numeric, \eqn{\tau} index of gene 1.}
#'   \item{\emph{class_dup1} Factor, expression-based class of gene 1.}
#'   \item{\emph{tau_dup2} Numeric, \eqn{\tau} index of gene 2.}
#'   \item{\emph{class_dup2} Factor, expression-based class of gene 2.}
#'   \item{\emph{tau_diff} Numeric, absolute \eqn{\tau} difference between
#'   duplicate pairs.}
#'   \item{\emph{class_preserved} Logical indicating if expression-based class
#'   is the same for both genes in a duplicate pair.}
#' }
#'
#' @examples
#' # Load data
#' data(ath_se)
#' data(ath_dups)
#' dups <- ath_dups
#' se <- ath_se
#' 
#' # Simulate different tissues
#' SummarizedExperiment::colData(se)$tissue <- c(rep("seed", 3), rep("leaf", 3))
#'
#' # Aggregate expression to tissue level
#' aexp <- aggregate_to_tissue(se, "tissue")
#' 
#' # Calculate Tau
#' tau_df <- calculate_tau(aexp)
#' 
#' # Compare Tau
#' ct <- compare_tau(tau_df, dups)
compare_tau <- function(tau_df, dups) {
    
    # Add Tau and class for each member of the pair
    tdf <- tau_df
    names(tdf)[2:3] <- c("tau_dup1", "class_dup1")
    fdups <- merge(dups, tdf, by = 1, all.x = TRUE)
    names(tdf)[2:3] <- c("tau_dup2", "class_dup2")
    fdups <- merge(fdups, tdf, by.x = 2, by.y = 1, all.x = TRUE)
    
    # Add columns `tau_diff` and `class_preserved`
    fdups$tau_diff <- abs(fdups$tau_dup1 - fdups$tau_dup2)
    fdups$class_preserved <- ifelse(
        fdups$class_dup1 == fdups$class_dup2, TRUE, FALSE
    )
    
    # Reorder columns to keep dup1 and dup2 in order
    cidx <- seq_len(ncol(fdups))
    cidx[c(1, 2)] <- cidx[c(2, 1)]
    fdups <- fdups[, cidx]
    
    return(fdups)
}



#' Compare coexpression module assignment for gene pairs
#'
#' @param dups A data frame of duplicate pairs, with duplicated genes
#' for each pair in columns 1 and 2. Additional columns are allowed, and will
#' be ignored if present.
#' @param genes_modules A 2-column data frame with genes in column 1
#' and coexpression modules to which they belong in column 1. 
#' This can be obtained, for instance, with \code{BioNERO::exp2gcn()}.
#'
#' @return A data frame with the same columns as the input to \strong{dups},
#' but with the following extra columns:
#' \itemize{
#'   \item{\emph{module_dup1} Character, name of module where gene 1 is.}
#'   \item{\emph{module_dup2} Character, name of module where gene 2 is.}
#'   \item{\emph{module_preservation} Factor, classification of the
#'   gene pair based on module preservation. Levels include 
#'   'preserved' (both genes in the same module), 'diverged' (genes
#'   in different modules), 'one_absent' (one of the genes is not any
#'   module), and 'both_absent' (both genes are not present in any module).}
#' }
#' 
#' @rdname compare_coex_modules
#' @export
#' @examples
#' data(ath_dups)
#' data(genes_and_modules)
#' dups <- ath_dups
#' genes_modules <- genes_and_modules
#' 
#' # Compare
#' cm <- compare_coex_modules(dups, genes_and_modules)
compare_coex_modules <- function(dups, genes_modules) {
    
    # Remove gray module and add module for each gene in a pair
    mod <- genes_modules
    names(mod)[2] <- "module_dup1"
    fdups <- merge(dups, mod, by = 1, all.x = TRUE)
    names(mod)[2] <- "module_dup2"
    fdups <- merge(fdups, mod, by.x = 2, by.y = 1, all.x = TRUE)
    fdups$module_dup1 <- gsub("grey|gray", NA, fdups$module_dup1)
    fdups$module_dup2 <- gsub("grey|gray", NA, fdups$module_dup2)
    
    # Add column `module_preserved`
    fdups$module_preservation <- apply(fdups, 1, function(x) {
        m <- c(x["module_dup1"], x["module_dup2"])
        nacount <- sum(is.na(m))
        
        if(nacount == 2) {
            p <- "both_absent"
        } else if(nacount == 1) {
            p <- "one_absent"
        } else if(m[1] == m[2]) {
            p <- "preserved"
        } else if(m[1] != m[2]) {
            p <- "diverged"
        }
    })
    fdups$module_preservation <- factor(
        fdups$module_preservation, 
        levels = c("preserved", "diverged", "one_absent", "both_absent")
    )
    
    # Reorder columns to keep dup1 and dup2 in order
    cidx <- seq_len(ncol(fdups))
    cidx[c(1, 2)] <- cidx[c(2, 1)]
    fdups <- fdups[, cidx]

    return(fdups)
}


#' Compare coexpression module eigengenes for diverged duplicate pairs
#'
#' @param module_comp A data frame with module assignment comparison
#' between duplicate pairs as returned by \code{compare_coex_modules()}.
#' @param ME A matrix or data frame of module eigengenes, with samples in rows
#' and eigengenes in columns. This can be obtained, for instance,
#' with \code{BioNERO::exp2gcn()}.
#' @param cor_method Character indicating the correlation method to use
#' to calculate pairwise correlations between module eigengenes.
#' Default: "spearman".
#'
#' @return A data frame with duplicate pairs (and any additional columns
#' with metadata on pairs) and a column named \emph{ME_cor} with correlation
#' between eigengenes of the modules where each gene is.
#' 
#' @details
#' This function aims to derive a continuous measure from a binary variable
#' (preserved vs diverged modules for genes in a duplicate pair). For that,
#' this function calculates correlation coefficients between module eigengenes
#' for diverged gene pairs. Such correlation coefficients can be interpreted
#' as a measure of how much genes in a pair are diverged (i.e., if they
#' are in very different modules, or in quite similar modules).
#' 
#' 
#' @importFrom stats cor
#' @export
#' @rdname compare_coex_me
#' @examples
#' data(ath_dups)
#' data(genes_and_modules)
#' data(ME)
#' dups <- ath_dups
#' genes_modules <- genes_and_modules
#' 
#' # Compare
#' module_comp <- compare_coex_modules(dups, genes_and_modules)
#' comp <- compare_coex_me(module_comp, ME)
compare_coex_me <- function(module_comp, ME, cor_method = "spearman") {
    
    # Remove grey module
    ME <- ME[, colnames(ME) != "MEgrey"]
    
    # Get pairwise correlations between module eigengenes
    me_cor <- cor(ME, method = cor_method)
    me_df <- cormat2long(me_cor, redundant = TRUE)
    
    # Get only pairs for which both genes are in modules
    fpairs <- module_comp[module_comp$module_preservation == "diverged", ]
    
    # Add column with correlation between MEs
    fpairs$m1m2 <- paste0(fpairs$module_dup1, "-", fpairs$module_dup2)
    cor_df <- data.frame(
        m1m2 = gsub("ME", "", paste0(me_df$var1, "-", me_df$var2)), 
        ME_cor = me_df$cor
    )
    fpairs <- merge(fpairs, cor_df, by = "m1m2")
    
    # Remove unnecessary columns
    remove <- c("m1m2", "module_dup1", "module_dup2", "module_preservation")
    fpairs[, remove] <- NULL
    
    return(fpairs)
}
