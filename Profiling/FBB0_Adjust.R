# Step 1: Extract OTU table (samples are rows, taxa are columns)
otu <- as(otu_table(ps.genus.S1), "matrix")

# Step 2: Extract FBB0 replicates
fbb0_1 <- otu["S1-FBB0-1", ]
fbb0_2 <- otu["S1-FBB0-2", ]

# Step 3: Extract sample metadata
sample_df <- as(sample_data(ps.genus.S1), "data.frame")

# Step 4: Normalize each sample by its matching FBB0 replicate with pseudocount
otu_norm <- otu
for (sample_id in rownames(otu)) {
  rep_val <- sample_df[sample_id, "Rep"]
  if (rep_val == 1) {
    otu_norm[sample_id, ] <- (otu[sample_id, ] + 0.01) / (fbb0_1 + 0.01)
  } else if (rep_val == 2) {
    otu_norm[sample_id, ] <- (otu[sample_id, ] + 0.01) / (fbb0_2 + 0.01)
  }
}

# Step 5: Remove taxa (columns) with any NaN, Inf, or -Inf
bad_taxa <- apply(otu_norm, 2, function(col) any(!is.finite(col)))
otu_norm_clean <- otu_norm[, !bad_taxa]

# Step 6: Remove FBB0 samples from both OTU matrix and sample data
fbb0_samples <- rownames(sample_df)[grepl("FBB0", rownames(sample_df))]
otu_final <- otu_norm_clean[!rownames(otu_norm_clean) %in% fbb0_samples, ]
sample_final <- sample_df[!rownames(sample_df) %in% fbb0_samples, ]

# Step 7: Rebuild phyloseq object with cleaned data
# Subset taxonomy table too
tax_final <- tax_table(ps.genus.S1)[colnames(otu_final), ]
otu_phyloseq <- otu_table(otu_final, taxa_are_rows = FALSE)
sample_phyloseq <- sample_data(sample_final)
ps.genus.S1.cleaned <- phyloseq(otu_phyloseq, tax_final, sample_phyloseq)

#Write the funciton
library(phyloseq)

normalize_by_fbb0 <- function(ps_object) {
  # Extract OTU table and metadata
  otu <- as(otu_table(ps_object), "matrix")
  if (taxa_are_rows(ps_object)) {
    otu <- t(otu)
  }
  
  sample_df <- as(sample_data(ps_object), "data.frame")
  
  # Identify FBB0 replicates
  fbb0_samples <- rownames(sample_df)[grepl("FBB0", rownames(sample_df))]
  fbb0_df <- sample_df[fbb0_samples, ]
  
  # Ensure there are exactly 2 replicates: Rep 1 and Rep 2
  if (!all(c(1, 2) %in% fbb0_df$Rep)) {
    stop("Both FBB0 Rep 1 and Rep 2 must be present.")
  }
  
  fbb0_1 <- otu[rownames(fbb0_df)[fbb0_df$Rep == 1], ]
  fbb0_2 <- otu[rownames(fbb0_df)[fbb0_df$Rep == 2], ]
  
  if (nrow(fbb0_1) != 1 || nrow(fbb0_2) != 1) {
    stop("Expected exactly one sample for each FBB0 replicate.")
  }
  
  # Normalize each sample by matching FBB0
  otu_norm <- otu
  for (sample_id in rownames(otu)) {
    rep_val <- sample_df[sample_id, "Rep"]
    if (rep_val == 1) {
      otu_norm[sample_id, ] <- (otu[sample_id, ] + 0.01) / (fbb0_1 + 0.01)
    } else if (rep_val == 2) {
      otu_norm[sample_id, ] <- (otu[sample_id, ] + 0.01) / (fbb0_2 + 0.01)
    }
  }
  
  # Remove taxa (columns) with any NaN/Inf/-Inf
  bad_taxa <- apply(otu_norm, 2, function(x) any(!is.finite(x)))
  otu_norm_clean <- otu_norm[, !bad_taxa]
  
  # Remove FBB0 samples
  keep_samples <- setdiff(rownames(otu_norm_clean), fbb0_samples)
  otu_final <- otu_norm_clean[keep_samples, ]
  sample_final <- sample_df[keep_samples, ]
  tax_final <- tax_table(ps_object)[colnames(otu_final), ]
  
  # Rebuild phyloseq object
  otu_phyloseq <- otu_table(otu_final, taxa_are_rows = FALSE)
  sample_phyloseq <- sample_data(sample_final)
  
  ps_cleaned <- phyloseq(otu_phyloseq, tax_final, sample_phyloseq)
  
  return(ps_cleaned)
}

