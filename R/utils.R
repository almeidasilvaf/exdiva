
#' Aggregate expression to tissue level
#'
#' @param se A `SummarizedExperiment` object with expression data in
#' TPM (or normalized by library size) and a `colData` slot.
#' @param tissue_column Character indicating the name of the column 
#' in \code{colData(se)} containing tissue information for each sample.
#' @param aggregation_method Character indicating the method to use to
#' aggregate samples from the same tissue before calculating \eqn{\tau}
#' (see details below). One of "median" or "mean". Default: "median".
#'
#' @return A gene expression matrix with genes in rows and tissues in columns.
#'
#' @importFrom SummarizedExperiment colData assay
#' @importFrom stats median
#' @rdname aggregate_to_tissue
#' @export
#' @examples
#' # Load data
#' data(ath_se)
#' se <- ath_se
#' 
#' # Simulate different tissues
#' SummarizedExperiment::colData(se)$tissue <- c(rep("seed", 3), rep("leaf", 3))
#'
#' # Aggregate expression
#' aexp <- aggregate_to_tissue(se, "tissue")
aggregate_to_tissue <- function(
        se, tissue_column, aggregation_method = "median"
) {
    
    # Aggregate expression per tissue
    cdata <- as.data.frame(colData(se))
    if(!tissue_column %in% names(cdata)) {
        stop("Column ", tissue_column, " was not found in the colData slot.")
    }
    tissues <- unique(cdata[[tissue_column]])
    
    exp_tissue <- Reduce(cbind, lapply(tissues, function(x) {
        
        ## Get expression for samples that belong to tissue {x}
        samples <- rownames(cdata[cdata[[tissue_column]] == x, ])
        exp_x <- assay(se)[, samples]
        
        ## Aggregate expression to tissue level
        exp_x <- as.matrix(apply(exp_x, 1, aggregation_method, na.rm = TRUE))
        
        return(exp_x)
    }))
    colnames(exp_tissue) <- tissues
    
    return(exp_tissue)
}

#' Classify genes into categories based on \eqn{\tau} indices
#'
#' @param aggregated_exp A gene expression matrix with genes in rows
#' and tissues in columns as returned by \code{aggregate_to_tissue}.
#' @param tau_df A 2-column data frame with variables \strong{gene}
#' and \strong{tau}.
#' 
#' @return A data frame as passed to \strong{tau_df}, but with an extra
#' column named `class`.
#' @noRd
classify_tau <- function(aggregated_exp, tau_df) {
    
    # Number of tissues in which each gene is expressed and stably expressed
    nexpressed <- apply(aggregated_exp, 1, function(x) sum(x >= 1))
    nstable <- apply(aggregated_exp, 1, function(x) sum(x >= 5))
    
    cdf <- cbind(tau_df, nexpressed, nstable)
    eclass <- apply(cdf, 1, function(x) {
        if(is.na(x["tau"])) {
            c <- NA
        } else if(x["nexpressed"] == 0) { # not expressed in any tissue
            c <- "Null"
        } else if(x["nstable"] == 0) { # not stably expressed in any tissue
            c <- "Weak"
        } else if(x["nstable"] > 0 & x["tau"] < 0.85) { # stable, Tau <0.85
            c <- "Broad"
        } else if(x["nstable"] > 0 & x["tau"] >= 0.85) { # stable, Tau >=0.85
            c <- "Specific"
        }
        
        return(c)
    })
    
    tau_df$class <- factor(
        eclass, levels = c("Null", "Weak", "Broad", "Specific")
    )
    
    return(tau_df)
}


#' Wrapper function to convert a symmetric matrix to a long data frame
#' 
#' @param cormat A correlation matrix (symmetric).
#' @param redundant Logical indicating whether to leave redundancy (upper
#' and lower triangle are the same).
#' 
#' @return A 3-column data frame with variables \strong{var1}, \strong{var2},
#' and \strong{cor}, representing the original correlation matrix in long
#' format.
#' 
#' @importFrom stats na.omit
#' @noRd
cormat2long <- function(cormat, redundant = FALSE) {
    
    # Get a data frame in long format (similar to pivot_longer())
    df <- cormat
    if(!redundant) {
        df[lower.tri(df, diag = TRUE)] <- NA
    }
    df <- na.omit(data.frame(as.table(df)))
    
    # Rename columns and fix data classes
    names(df) <- c("var1", "var2", "cor")
    df$var1 <- as.character(df$var1)
    df$var2 <- as.character(df$var2)
    df$cor <- as.numeric(df$cor)
    
    return(df)
}


