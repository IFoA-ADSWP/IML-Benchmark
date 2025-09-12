to_num <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", "", x)               # strip thousand separators if any
  suppressWarnings(as.numeric(x))     # coerce; NAs where non-numeric
}

train_GAM <- function(dt){
  fac_cols <- c("Area","VehBrand","VehGas","Region")
  num_cols <- c("VehPower","VehAge","DrivAge","BonusMalus","Density")
  for (v in intersect(fac_cols, names(dt))) dt[[v]] <- as.factor(dt[[v]])
  for (v in intersect(num_cols, names(dt))) dt[[v]] <- to_num(dt[[v]])
  
  gam_formula <- ClaimNb ~ s(DrivAge, k = 10) 
  + s(VehAge, k = 10)
  + s(VehPower, k = 8) 
  + s(BonusMalus, k = 10) 
  + s(Density, k = 10) 
  + te(DrivAge,BonusMalus, k = c(8, 8)) 
  + s(Area, bs = "re") 
  + s(VehBrand, bs = "re") 
  + s(VehGas, bs = "re") 
  + s(Region, bs = "re")
  
  mgcv::bam(
    formula  = gam_formula,
    family   = poisson(link="log"),
    data     = dt,
    method   = "fREML",
    discrete = TRUE
    # , select = TRUE   # optional: shrinkage/feature selection
  )
}

predict_GAM <- function(model, newdata, type="response"){
  fac_cols <- c("Area","VehBrand","VehGas","Region")
  num_cols <- c("VehPower","VehAge","DrivAge","BonusMalus","Density")
  
  # numeric safely
  for (v in intersect(num_cols, names(newdata))) newdata[[v]] <- to_num(newdata[[v]])
  
  # align factor levels to training levels; unseen -> NA (or map to "Other" if available)
  for (v in intersect(fac_cols, names(newdata))) {
    trained_lvls <- model$xlevels[[v]]
    if (is.null(trained_lvls)) trained_lvls <- levels(model$model[[v]])
    newdata[[v]] <- factor(newdata[[v]], levels = trained_lvls)
    # optional: if you trained with an "Other" level, map NAs to it:
    # if ("Other" %in% trained_lvls) newdata[[v]][is.na(newdata[[v]])] <- "Other"
  }
  
  as.vector(predict(model, newdata=newdata, type=type))
}
