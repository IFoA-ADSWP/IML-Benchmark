source("init.R")
source("py_init.R")
source("Models/train_GAM.R")
source("Models/train_EBM.R")
source("Models/train_localglmnet.R")

n_test_subset <- 100000
plot_var <- "DrivAge"
out_dir <- "Results/IML_charts_CV5"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
idx_sub <- sample(seq_len(nrow(test_df)), size = n_test_subset, replace = FALSE)

test_df_plot <- test_df[idx_sub, , drop = FALSE]
x_test <- dplyr::select(test_df_plot, -ClaimNb)

# GLM PDP via iml -------------------------------------------------------------
glm_pred <- Predictor$new(
  model = models$CV_5$glm_model,
  data = x_test,
  y = test_df_plot$ClaimNb,
  predict.function = function(model, newdata) {
    as.numeric(stats::predict(model, newdata = newdata, type = "response"))
  }
)

p_glm <- plot(FeatureEffect$new(glm_pred, feature = plot_var, method = "pdp")) +
  labs(title = "GLM PDP - DrivAge", x = "DrivAge", y = "Predicted claim frequency") +
  gg_style

ggsave(file.path(out_dir, "glm_pdp_drivage.png"), p_glm, width = 8, height = 5, dpi = 140)

# IBLM charts from package docs ------------------------------------------------
iblm_ex <- IBLM::explain_iblm(models$CV_5$IBLM, test_df_plot)

p_iblm_scatter <- iblm_ex$beta_corrected_scatter(plot_var) +
  labs(title = "IBLM beta-corrected scatter - DrivAge") +
  gg_style

ggsave(file.path(out_dir, "iblm_beta_scatter_drivage.png"), p_iblm_scatter, width = 8, height = 5, dpi = 140)

p_iblm_density <- iblm_ex$beta_corrected_density(plot_var) +
  labs(title = "IBLM beta-corrected density - DrivAge") +
  gg_style

ggsave(file.path(out_dir, "iblm_beta_density_drivage.png"), p_iblm_density, width = 8, height = 5, dpi = 140)

# LocalGLMnet charts (moved helpers into train_localglmnet.R) -----------------
x_local <- transform_localglmnet_input(x_test, models$CV_5$LocalGLMnet$encoder)

p_local_att <- plot_localglmnet_attention(
  model = models$CV_5$LocalGLMnet,
  X = x_local,
  feature_names = colnames(x_local),
  plot_vars = plot_var,
  encoder = models$CV_5$LocalGLMnet$encoder,
  raw_df = x_test,
  sample_size = nrow(x_local),
  ncol = 1
) +
  labs(title = "LocalGLMnet attention - DrivAge") +
  gg_style

ggsave(file.path(out_dir, "localglmnet_attention_drivage.png"), p_local_att, width = 8, height = 5, dpi = 140)

p_local_contrib <- plot_localglmnet_contributions(
  model = models$CV_5$LocalGLMnet,
  X = x_local,
  feature_names = colnames(x_local),
  plot_vars = plot_var,
  encoder = models$CV_5$LocalGLMnet$encoder,
  raw_df = x_test,
  sample_size = nrow(x_local),
  ncol = 1
) +
  labs(title = "LocalGLMnet contribution - DrivAge") +
  gg_style

ggsave(file.path(out_dir, "localglmnet_contribution_drivage.png"), p_local_contrib, width = 8, height = 5, dpi = 140)

# GAM chart in ggplot via iml PDP ---------------------------------------------
gam_pred <- Predictor$new(
  model = models$CV_5$GAM_model,
  data = x_test,
  y = test_df_plot$ClaimNb,
  predict.function = function(model, newdata) {
    predict_GAM(model, newdata, type = "response")
  }
)

p_gam <- plot(FeatureEffect$new(gam_pred, feature = plot_var, method = "pdp")) +
  labs(title = "GAM PDP - DrivAge", x = "DrivAge", y = "Predicted claim frequency") +
  gg_style

ggsave(file.path(out_dir, "gam_pdp_drivage.png"), p_gam, width = 8, height = 5, dpi = 140)

# EBM one-way style chart in R -------------------------------------------------
# We can directly extract per-term additive contributions with eval_terms.
term_names <- py_to_r(models$CV_5$EBM$model$term_names_)
term_idx <- match(plot_var, term_names)

if (!is.na(term_idx)) {
  ebm_terms <- py_to_r(models$CV_5$EBM$model$eval_terms(r_to_py(data.matrix(x_test))))
  ebm_df <- data.frame(
    DrivAge = x_test[[plot_var]],
    effect = ebm_terms[, term_idx]
  )

  p_ebm <- ggplot(ebm_df, aes(x = DrivAge, y = effect)) +
    geom_point(alpha = 0.18, color = "#2D6A4F") +
    geom_smooth(method = "loess", se = TRUE, color = "#1B4332", fill = "#95D5B2") +
    labs(title = "EBM one-way effect proxy - DrivAge", x = "DrivAge", y = "Additive term effect") +
    gg_style

  ggsave(file.path(out_dir, "ebm_oneway_drivage.png"), p_ebm, width = 8, height = 5, dpi = 140)
} else {
  message("EBM term for DrivAge not found; skipping one-way plot.")
}

# XGBoost SHAP via predict(..., predcontrib=TRUE) + shapviz -------------------
xgb_shap <- predict(
  models$CV_5$XGB,
  xgboost::xgb.DMatrix(data.matrix(x_test)),
  predcontrib = TRUE
)

xgb_shap <- xgb_shap[, colnames(xgb_shap) != "(Intercept)", drop = FALSE]
sv <- shapviz::shapviz(xgb_shap, X = as.data.frame(x_test))

p_xgb <- shapviz::sv_dependence(sv, v = plot_var, color_var = NULL) +
  labs(title = "XGBoost SHAP dependence - DrivAge", x = "DrivAge", y = "SHAP value") +
  gg_style

ggsave(file.path(out_dir, "xgb_shap_drivage.png"), p_xgb, width = 8, height = 5, dpi = 140)

message("Saved charts to: ", out_dir)
