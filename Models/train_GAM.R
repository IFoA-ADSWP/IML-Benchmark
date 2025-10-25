to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", "", x)               # strip thousand separators if any
  suppressWarnings(as.numeric(x))     # coerce; NAs where non-numeric
}

train_GAM <- function(dt){
  fac_cols <- c("Area", "VehBrand", "VehGas", "Region")
  num_cols <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")
  
  # Convert variables to correct types
  for (v in intersect(fac_cols, names(dt))) dt[[v]] <- as.factor(dt[[v]])
  for (v in intersect(num_cols, names(dt))) dt[[v]] <- to_num(dt[[v]])
  
  # GAM formula (frequency model with numeric-only interactions)
  gam_formula <- ClaimNb ~ 
    s(DrivAge, k = 10) +
    s(VehAge, k = 10) +
    s(VehPower, k = 10) +
    s(BonusMalus, k = 10) +
    s(Density, k = 10) +
    s(Area, bs = "re") +
    s(VehBrand, bs = "re") +
    s(VehGas, bs = "re") +
    s(Region, bs = "re")
  
  # Fit GAM model
  mgcv::bam(
    formula  = gam_formula,
    family   = poisson(link = "log"),
    data     = dt,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE     # auto shrinkage / variable selection
  )
}

predict_GAM <- function(model, newdata, type="response"){
  fac_cols <- c("Area", "VehBrand", "VehGas", "Region")
  num_cols <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")
  
  # Convert numeric safely
  for (v in intersect(num_cols, names(newdata))) newdata[[v]] <- to_num(newdata[[v]])
  
  # Align factor levels to training data
  for (v in intersect(fac_cols, names(newdata))) {
    trained_lvls <- model$xlevels[[v]]
    if (is.null(trained_lvls)) trained_lvls <- levels(model$model[[v]])
    newdata[[v]] <- factor(newdata[[v]], levels = trained_lvls)
  }
  
  as.vector(predict(model, newdata = newdata, type = type))
}
