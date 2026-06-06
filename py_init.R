library(reticulate)
reticulate::use_virtualenv(paste0(getwd(),"/.venv"), required = TRUE)
# interpret = reticulate::import("interpret.glassbox")


if(!py_module_available("tensorflow")){
    tensorflow::install_tensorflow(method = "virtualenv", #
    version = "default", 
    envname = paste0(getwd(),"/.venv"))
}
keras::is_keras_available()
seed = 2

tensorflow::set_random_seed(seed, disable_gpu = FALSE)
py_config()
tensorflow::tf_version()


