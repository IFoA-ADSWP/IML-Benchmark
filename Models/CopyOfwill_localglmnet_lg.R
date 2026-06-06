# Sys.setenv(TF_ENABLE_ONEDNN_OPTS = "0")
# library(reticulate)
# library(keras3)
source("py_init.R")
source("init.R")

freMTPL2freq <- dt_list$fre_mtpl2_freq 
numeric_cols <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")
X <- as.matrix(freMTPL2freq[, numeric_cols])
Y <- freMTPL2freq$ClaimNb
complete_idx <- complete.cases(X, Y)
X <- X[complete_idx, ]
Y <- Y[complete_idx]
n <- nrow(X)/2

X_L <- X[1:n, ];         Y_L <- Y[1:n]
X_T <- X[(n+1):(2*n), ]; Y_T <- Y[(n+1):(2*n)]

cat("Learning set:", nrow(X_L), "rows\n")
cat("Test set:    ", nrow(X_T), "rows\n")
cat("Columns in X_L:", ncol(X_L), "\n")  # creating the train/test sets slightly differently

# ── Build LocalGLMnet ──────────────────────────────────────────────────────────
build_localglmnet <- function(q = 5) {
  
  Design <- layer_input(shape = c(q), dtype = "float32", name = "Design")
  
  Attention <- Design %>%
    layer_dense(units = 20, activation = "tanh") %>%
    layer_dense(units = 15, activation = "tanh") %>%
    layer_dense(units = 10, activation = "tanh") %>%
    layer_dense(units = q,  activation = "linear", name = "Attention")
  
  # Element-wise product then sum via Dense with fixed weights
  Response <- layer_multiply(list(Design, Attention)) %>%
    layer_dense(units = 1, activation = "linear", name = "Response")
  
  model <- keras_model(inputs = Design, outputs = Response)
  
  model %>% compile(loss = "poisson", 
                    optimizer = optimizer_adam())
  model
}

model <- build_localglmnet(q = 5)   # <-- must match q above
summary(model)

# ── Fit ────────────────────────────────────────────────────────────────────────
history <- model %>% fit(
  x                = X_L,
  y                = Y_L,
  epochs           = 100,
  batch_size       = 50000,
  validation_split = 0.2,
  callbacks        = list(
    callback_early_stopping(
      monitor              = "val_loss",
      patience             = 10,
      restore_best_weights = TRUE
    )
  ),
  verbose = 1
)


paste0("In-sample  MSE:", model %>%  evaluate(X_L, Y_L, verbose = 0), "\n")
paste0("Out-sample MSE:", model %>% evaluate(X_T, Y_T, verbose = 0), "\n")

library(ggplot2)
library(tidyverse)

# ── Extract attention weights ──────────────────────────────────────────────────
attention_layer <- keras_model(inputs = model$input, 
                               outputs = model$get_layer("Attention")$output)
attention_weights <- predict(attention_layer, X_T)

# ── Prepare data for plotting ──────────────────────────────────────────────────
set.seed(42)
sample_idx <- sample(1:nrow(X_T), size = min(5000, nrow(X_T)), replace = FALSE)

plot_data <- data.frame(
  feature_1 = X_T[sample_idx, 1],
  feature_2 = X_T[sample_idx, 2],
  feature_3 = X_T[sample_idx, 3],
  feature_4 = X_T[sample_idx, 4],
  feature_5 = X_T[sample_idx, 5],
  attention_1 = attention_weights[sample_idx, 1],
  attention_2 = attention_weights[sample_idx, 2],
  attention_3 = attention_weights[sample_idx, 3],
  attention_4 = attention_weights[sample_idx, 4],
  attention_5 = attention_weights[sample_idx, 5]
)

# ── Feature names ──────────────────────────────────────────────────────────────
feature_names <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")

# ── Convert to long format ─────────────────────────────────────────────────────
plot_data_long <- plot_data %>%
  pivot_longer(cols = starts_with("feature_"), 
               names_to = "feature", names_prefix = "feature_", 
               values_to = "feature_value") %>%
  pivot_longer(cols = starts_with("attention_"), 
               names_to = "attention_feature", names_prefix = "attention_",
               values_to = "attention_value") %>%
  filter(feature == attention_feature) %>%
  mutate(feature_num = as.numeric(feature),
         feature_name = factor(feature_num, 
                               levels = 1:5, 
                               labels = feature_names))

# ── Create faceted plot with feature names ─────────────────────────────────────
combined_plot <- ggplot(plot_data_long, aes(x = feature_value, y = attention_value)) +
  geom_point(alpha = 0.3, size = 1) +
  facet_wrap(~feature_name, scales = "free", ncol = 3, 
             labeller = labeller(feature_name = c(VehPower = "VehPower (1)",
                                                  VehAge = "VehAge (2)",
                                                  DrivAge = "DrivAge (3)",
                                                  BonusMalus = "BonusMalus (4)",
                                                  Density = "Density (5)"))) +
  labs(title = "Regression Attentions",
       x = "Feature values",
       y = "Regression attention β̂ⱼ(x)") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        strip.text = element_text(size = 10, face = "bold"))

print(combined_plot)


library(ggplot2)
library(tidyverse)

# ── Extract attention weights ──────────────────────────────────────────────────
attention_layer <- keras_model(inputs = model$input, 
                               outputs = model$get_layer("Attention")$output)
attention_weights <- predict(attention_layer, X_T)

# ── Prepare data for plotting ──────────────────────────────────────────────────
set.seed(42)
sample_idx <- sample(1:nrow(X_T), size = min(5000, nrow(X_T)), replace = FALSE)

# Calculate feature contributions: attention_weight * feature_value
contributions <- data.frame(
  feature_1 = attention_weights[sample_idx, 1] * X_T[sample_idx, 1],
  feature_2 = attention_weights[sample_idx, 2] * X_T[sample_idx, 2],
  feature_3 = attention_weights[sample_idx, 3] * X_T[sample_idx, 3],
  feature_4 = attention_weights[sample_idx, 4] * X_T[sample_idx, 4],
  feature_5 = attention_weights[sample_idx, 5] * X_T[sample_idx, 5],
  x_1 = X_T[sample_idx, 1],
  x_2 = X_T[sample_idx, 2],
  x_3 = X_T[sample_idx, 3],
  x_4 = X_T[sample_idx, 4],
  x_5 = X_T[sample_idx, 5]
)

# ── Feature names ──────────────────────────────────────────────────────────────
feature_names <- c("VehPower", "VehAge", "DrivAge", "BonusMalus", "Density")

# ── Convert to long format ─────────────────────────────────────────────────────
plot_data_long <- contributions %>%
  pivot_longer(cols = starts_with("feature_"), 
               names_to = "feature", names_prefix = "feature_", 
               values_to = "contribution") %>%
  pivot_longer(cols = starts_with("x_"), 
               names_to = "x_feature", names_prefix = "x_",
               values_to = "feature_value") %>%
  filter(feature == x_feature) %>%
  mutate(feature_num = as.numeric(feature),
         feature_name = factor(feature_num, 
                               levels = 1:5, 
                               labels = feature_names))

# ── Create faceted plot ────────────────────────────────────────────────────────
combined_plot <- ggplot(plot_data_long, aes(x = feature_value, y = contribution)) +
  geom_point(alpha = 0.3, size = 1, color = "black") +
  geom_smooth(method = "loess", se = FALSE, color = "red", linewidth = 1) +
  facet_wrap(~feature_name, scales = "free", ncol = 3,
             labeller = labeller(feature_name = c(VehPower = "feature contribution: feature x1",
                                                  VehAge = "feature contribution: feature x2",
                                                  DrivAge = "feature contribution: feature x3",
                                                  BonusMalus = "feature contribution: feature x4",
                                                  Density = "feature contribution: feature x5"))) +
  # Add reference lines (mean, quartiles)
  geom_hline(yintercept = 0, color = "orange", linewidth = 1, linetype = "solid") +
  labs(title = "Feature Contributions",
       x = "Feature values xⱼ",
       y = "Contribution β̂ⱼ(x) × xⱼ") +
  theme_minimal() +
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        strip.text = element_text(size = 9, face = "bold"))

print(combined_plot)

ggsave("feature_contributions.png", combined_plot, width = 14, height = 8, dpi = 300)