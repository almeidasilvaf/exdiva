Data acquisition
================

# Data in data/

Here, we will describe the code used to create the example data sets in
the `data/` directory of this package.

``` r
set.seed(123) # for reproducibility
library(tidyverse)
library(BioNERO)
```

## `ath_dups.rda`

This file contains a data frame with duplicate gene pairs in the genome
of *Arabidopsis thaliana*, and it was obtained from doubletroubledb.

``` r
# Load download .rda file
load("~/Downloads/plants_duplicates.rda")

# Get data frame of duplicates for A. thaliana and save it
ath_dups <- plants_duplicates$arabidopsis_thaliana |>
    dplyr::mutate(
        dup1 = stringr::str_replace_all(dup1, "ara_", ""),
        dup2 = stringr::str_replace_all(dup2, "ara_", "")
    )

rownames(ath_dups) <- NULL

usethis::use_data(ath_dups, compress = "xz")
```

## `ath_se.rda`

This file contains a `SummarizedExperiment` object with RNA-seq data for
*Arabidopsis thaliana*, publicly available under SRA Project SRP213876.
Expression data and sample metadata were obtained from Refine Bio.

``` r
data("ath_dups")

# Get expression matrix (in TPM)
exp <- readr::read_tsv(
    "~/Downloads/SRP213876_refinebio/SRP213876/SRP213876.tsv",
    show_col_types = FALSE
) |>
    column_to_rownames("Gene") |>
    as.matrix()

# Filter out singletons
keep <- intersect(rownames(exp), unique(c(ath_dups$dup1, ath_dups$dup2)))
exp <- exp[keep, ]

# Get sample metadata
coldata <- read_tsv(
    "~/Downloads/SRP213876_refinebio/SRP213876/metadata_SRP213876.tsv",
    show_col_types = FALSE
) |>
    mutate(
        treatment = case_when(
            str_detect(refinebio_title, "FRp") ~ "Far-red",
            TRUE ~ "Red"
        )
    ) |>
    select(
        experiment = experiment_accession,
        accession = refinebio_accession_code,
        tissue = refinebio_subject,
        treatment
    ) |>
    tibble::column_to_rownames("accession")

# Create SummarizedExperiment object
ath_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(TPM = exp),
    colData = coldata
)

# Save object
usethis::use_data(ath_se, compress = "xz", overwrite = TRUE)
```

## `genes_and_modules.rda` and `ME.rda`

To create these objects, we inferred a gene coexpression network with
*[BioNERO](https://bioconductor.org/packages/3.18/BioNERO)* using data
from SRP201971. Expression data and metadata were obtained from Refine
Bio.

``` r
data(ath_dups)

exp <- readr::read_tsv(
    "~/Downloads/SRP201971_refinebio/SRP201971/SRP201971.tsv",
    show_col_types = FALSE
) |>
    column_to_rownames("Gene") |>
    as.matrix()

# Filter out singletons
keep <- intersect(rownames(exp), unique(c(ath_dups$dup1, ath_dups$dup2)))
exp <- exp[keep, ]

# Get sample metadata
coldata <- read_tsv(
    "~/Downloads/SRP201971_refinebio/SRP201971/metadata_SRP201971.tsv",
    show_col_types = FALSE
) |>
    mutate(
        treatment = case_when(
            str_detect(refinebio_title, "H2O") ~ "H2O",
            str_detect(refinebio_title, "flg22") ~ "flg22",
            str_detect(refinebio_title, "nlp20") ~ "nlp20",
            str_detect(refinebio_title, "C6") ~ "C6"
        )
    ) |>
    select(
        experiment = experiment_accession,
        accession = refinebio_accession_code,
        tissue = refinebio_subject,
        time = refinebio_time,
        treatment
    ) |>
    tibble::column_to_rownames("accession")

# Create SummarizedExperiment object
se_mamp <- SummarizedExperiment::SummarizedExperiment(
    assays = list(TPM = exp),
    colData = coldata
)

# Filter:
# 1) Keep only genes with TPM >= 1 in at least 20% of the samples
# 2) PC-based correction for confounders
filt_se <- exp_preprocess(
    se_mamp, 
    min_exp = 1,
    method = "percentage",
    min_percentage_samples = 0.2,
    Zk_filtering = FALSE, 
    cor_method = "pearson"
)

# Apply scale-free topology fit
WGCNA::allowWGCNAThreads(nThreads = 4)
sft <- SFT_fit(filt_se, net_type = "signed", cor_method = "pearson")

# Infer GCN
gcn <- exp2gcn(
    filt_se, 
    net_type = "signed", 
    SFTpower = sft$power,
    cor_method = "pearson"
)

# Get objects
ME <- gcn$MEs
genes_and_modules <- gcn$genes_and_modules

usethis::use_data(ME, compress = "xz")
usethis::use_data(genes_and_modules, compress = "xz")
```
