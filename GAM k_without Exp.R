source("init.R")
source("Models/train_GAM.R")

dt_list$fre_mtpl2_freq$Exposure=1

CV <- 5
set.seed(1)

CV_vec <- sample(1:CV, replace = TRUE, size = nrow(dt_list$fre_mtpl2_freq))

models  <- list()
results <- list()

# bring homog back
losses <- data.frame(
  CV    = paste0("CV_", 1:CV),
  homog = NA_real_,
  glm   = NA_real_,
  GAM   = NA_real_,
  GAM_a = NA_real_,
  GAM_b = NA_real_,
  GAM_c = NA_real_
)

for (i in 1:CV) {
  train_rows <- which(CV_vec != i)
  iter       <- paste0("CV_", i)
  
  models[[iter]] <- list()
  
  train_df <- dt_list$fre_mtpl2_freq[train_rows, ]
  test_df  <- dt_list$fre_mtpl2_freq[-train_rows, ]
  
  # drop only IDpol
  train_df$IDpol <- NULL
  test_df$IDpol  <- NULL
  
  # test df
  results[[iter]] <- data.frame(
    ID     = dt_list$fre_mtpl2_freq$IDpol[-train_rows],
    actual = dt_list$fre_mtpl2_freq$ClaimNb[-train_rows]
  )
  
  ## ---------- homog ----------
  lambda_bar <- sum(train_df$ClaimNb) / sum(train_df$Exposure)   # global rate
  results[[iter]]$homog <- lambda_bar * test_df$Exposure
  losses$homog[i] <- poisson_deviance(
    y_true = results[[iter]]$actual,
    y_pred = results[[iter]]$homog
  )
  
  ## ---------- GLM ----------
  models[[iter]]$glm_model <- glm(
    ClaimNb ~DrivAge+VehAge+BonusMalus+VehPower +Density +Area + VehBrand+VehGas+Region,
    family = poisson,
    offset = log(Exposure),
    data   = train_df
  )
  
  results[[iter]]$glm <- predict(
    models[[iter]]$glm_model,
    test_df,
    type = "response"
  )
  
  losses$glm[i] <- poisson_deviance(
    y_true = results[[iter]]$actual,
    y_pred = results[[iter]]$glm
  )
  
  ## ---------- GAM ----------
  models[[iter]]$GAM_model   <- train_GAM(train_df)
  results[[iter]]$GAM        <- predict_GAM(models[[iter]]$GAM_model,   test_df, type = "response")
  losses$GAM[i] <- poisson_deviance(results[[iter]]$actual, results[[iter]]$GAM)
  
  ## ---------- GAM_a ----------
  models[[iter]]$GAM_a_model <- train_GAM_a(train_df)
  results[[iter]]$GAM_a      <- predict_GAM(models[[iter]]$GAM_a_model, test_df, type = "response")
  losses$GAM_a[i] <- poisson_deviance(results[[iter]]$actual, results[[iter]]$GAM_a)
  
  ## ---------- GAM_b ----------
  models[[iter]]$GAM_b_model <- train_GAM_b(train_df)
  results[[iter]]$GAM_b      <- predict_GAM(models[[iter]]$GAM_b_model, test_df, type = "response")
  losses$GAM_b[i] <- poisson_deviance(results[[iter]]$actual, results[[iter]]$GAM_b)
  
  ## ---------- GAM_c ----------
  models[[iter]]$GAM_c_model <- train_GAM_c(train_df)
  results[[iter]]$GAM_c      <- predict_GAM(models[[iter]]$GAM_c_model, test_df, type = "response")
  losses$GAM_c[i] <- poisson_deviance(results[[iter]]$actual, results[[iter]]$GAM_c)
}

# overall and per-fold results (your original style)
poiss_per_CV <- rbind(
  losses,
  losses %>%
    tidyr::pivot_longer(cols = !CV) %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(mean_poiss = mean(value), .groups = "drop") %>%
    dplyr::arrange(mean_poiss) %>%
    tidyr::pivot_wider(values_from = mean_poiss, names_from = name) %>%
    dplyr::mutate(CV = "mean_poiss")
)

# pinball
poiss_per_CV %>%
  dplyr::mutate_if(
    is.numeric,
    ~ dplyr::if_else(. == homog, ., 1 - . / homog)
  ) %>%
  dplyr::select(-homog) %>%
  dplyr::mutate_if(
    is.numeric,
    ~ scales::percent(., accuracy = 0.1)
  )

library('mgcv')
gam.check(models$CV_1$GAM_model)
gam.check(models$CV_1$GAM_a_model)
gam.check(models$CV_1$GAM_b_model)
gam.check(models$CV_1$GAM_c_model)


