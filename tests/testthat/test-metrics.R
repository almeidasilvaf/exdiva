
# Load data ----
data(ath_se)
data(ath_dups)
data(genes_and_modules)
data(ME)

dups <- ath_dups[ath_dups$type == "TD", ]
se <- ath_se

# Start tests ----
test_that("calculate_cor() calculates correlation between gene pairs", {
    
    cor1 <- calculate_cor(se, dups[1:5, ], method = "pearson")
    
    expect_equal(class(cor1), "data.frame")
    expect_true("cor" %in% names(cor1))
})

test_that("calculate_tau() and compare_tau() returns data frames", {
    
    se2 <- se
    SummarizedExperiment::colData(se2)$tissue <- c(
        rep("seed", 3), rep("leaf", 3)
    )
    aexp <- aggregate_to_tissue(se2, "tissue")
    tau_df <- calculate_tau(aexp)
    
    ctau <- compare_tau(tau_df, dups)
    
    expect_error(aggregate_to_tissue(se2, "nonexistent_columm"))
    
    expect_equal(class(tau_df), "data.frame")
    expect_true("tau" %in% names(tau_df))
    expect_true("class" %in% names(tau_df))
    
    expect_equal(class(ctau), "data.frame")
    expect_true("class_preserved" %in% names(ctau))
})

test_that("compare_coex_modules() and compare_coex_me() return data frames", {
    
    cm <- compare_coex_modules(dups, genes_and_modules)
    cme <- compare_coex_me(cm, ME)
    
    r <- cormat2long(cor(ME))
    
    expect_equal(class(cm), "data.frame")
    expect_equal(class(cme), "data.frame")
    expect_true("module_preservation" %in% names(cm))
    expect_true("ME_cor" %in% names(cme))
})
