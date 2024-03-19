
# Load data ----
data(ath_se)

# Start tests ----
test_that("cor_nulls() returns a vector of nulls", {
    n <- cor_nulls(ath_se, subset = rownames(ath_se), n = 10)
    
    expect_equal(length(n), 10)
    expect_equal(class(n), "numeric")
})
