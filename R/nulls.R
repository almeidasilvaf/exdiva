
#' Generate a null distribution of correlation coefficients
#' 
#' @param se A `SummarizedExperiment` object with expression data.
#' @param subset (Optional) Character vector with genes to subset from the
#' expression matrix.
#' @param method Character indicating which correlation coefficient to compute.
#' One of 'pearson', 'spearman', or 'kendall'. Default: 'spearman'.
#' @param n Numeric indicating how many pairs to randomly sample. Default: 1e4.
#'
#' @return A numeric vector with a null distribution.
#' 
#' @noRd
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor
cor_nulls <- function(se, subset = NULL, method = "spearman", n = 1e4) {
    
    exp <- assay(se)
    genes <- rownames(exp)
    
    # Remove genes with zero standard deviation
    zero_sd <- apply(exp, 1, sd)
    zero_sd <- names(zero_sd[zero_sd == 0])
    if(length(zero_sd) > 0) { 
        genes <- genes[!genes %in% zero_sd]
    }
    if(!is.null(subset)) { genes <- genes[genes %in% subset] }
    
    # Sample `n` pairs with replacement
    nulls <- unlist(lapply(seq_len(n), function(x) {
        
        sg <- sample(genes, 2, replace = FALSE)
        corr <- cor(exp[sg[1], ], exp[sg[2], ], method = method)
        
        return(corr)
    }))
    
    return(nulls)
}

