train_df=read.csv("freMTPL2freq.csv")
# --- Data preprocessing before fitting GAM ---
fac_cols <- c("Area", "VehBrand", "VehGas", "Region")
num_cols <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")

# Convert to correct types
for (v in intersect(fac_cols, names(train_df))) {
  train_df[[v]] <- as.character(train_df[[v]])
  train_df[[v]][is.na(train_df[[v]]) | train_df[[v]] == ""] <- "Missing"
  train_df[[v]] <- as.factor(train_df[[v]])
}

for (v in intersect(num_cols, names(train_df))) {
  # safely convert to numeric
  train_df[[v]] <- suppressWarnings(as.numeric(as.character(train_df[[v]])))
  # simple imputation for missing values
  train_df[[v]][is.na(train_df[[v]])] <- median(train_df[[v]], na.rm = TRUE)
}

# --- Fit GAM model (simple & robust) ---
gam_fit <- mgcv::gam(
  ClaimNb ~ 
    s(DrivAge, k = 10) +
    s(VehAge, k = 10) +
    s(VehPower, k = 10) +
    s(BonusMalus, k = 10) +
    s(Density, k = 10) +
    s(Area, bs = "re") +
    s(VehBrand, bs = "re") +
    s(VehGas, bs = "re") +
    s(Region, bs = "re"),
  data   = train_df,
  family = poisson(link = "log"),
  method = "REML",
  select = TRUE    # automatic feature shrinkage
)
