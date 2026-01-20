library(reticulate)
python_path <- "/usr/bin/python3"
venv_path <- paste0(getwd(),"/venv")
if (!dir.exists(venv_path)) {dir.create(venv_path)}
reticulate::virtualenv_create(venv_path, python = python_path)
reticulate::use_virtualenv(venv_path, required = TRUE)

# Install python packages
reticulate::py_install("interpret", envname = venv_path)
reticulate::py_install("optuna", envname = venv_path)

interpret <- reticulate::import("interpret.glassbox")
optuna <- reticulate::import("optuna")
