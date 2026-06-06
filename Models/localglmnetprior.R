



train_LocalGLMnet = function(dt,
                             y,
                             vdt, #validation set
                             n_layers = 3, #Layers excluding the input layer, the output layer and the GLM.
                             n_units = c(20, 15, 10), #dim of the Layers excluding the input layer, the output layer and the GLM.
                             n_activ = c('tanh','tanh','tanh'),
                             n_dropout = c(F,F),
                             lr = 0.001, #learning rate
                             bs = 2^14, #batch size
                             ep = 500, #epochs
                             verbo = 0 #verbose
){
  
  log_size_hom <- log(mean(y))
  
  if(var(c(n_layers,length(n_units),length(n_activ))) != 0){stop("inconsistent n_h_layers, n_activ, n_units")}
  
  # neural network structure
  Design  <- layer_input(shape = c(dim(dt)[2]), dtype = 'float32', name = 'design') #Input layer
  
  # "Hidden layers" - all the layer excluding the input layer, the output layer and the glm
  for (i in 1:length(n_units)) {
    
    if (i==1) {
      Attention <- Design %>%    
        layer_dense(units = n_units[i], activation = n_activ[i], name=paste0('layer', (i)))
    }
    else if (i>1) {
      Attention = Attention %>%    
        layer_dense(units = n_units[i], activation = n_activ[i], name=paste0('layer', (i)))
    }
    
    if(n_dropout[i]){Attention = Attention %>% layer_dropout(rate = 0.2)}
    
  }
  
  #output layer definition
  Attention =Attention %>%
    layer_dense(units=c(dim(dt)[2]), activation='linear', name='attention')
  
  #GLM with fixed loglink function
  Output <- list(Design, Attention) %>% layer_dot(name='LocalGLM', axes=1) %>% 
    layer_dense(
      units=1, activation='exponential', name='output',
      weights=list(array(0, dim=c(1,1)), array(log_size_hom, dim=c(1)))
    )
  
  model_localGLMnet <- keras_model(inputs = list(Design), outputs = c(Output))

  model_localGLMnet %>% 
    compile(
    loss = custom_poisson,
    optimizer= optimizer_adam(learning_rate = lr)
  ) 
  
  fit_LocalGLMnet <- model_localGLMnet %>% 
    fit(
      list(as.matrix(dt)), list(as.matrix(as.numeric(y))),
      validation_data = vdt,
      batch_size = bs,
      epochs = ep,
      shuffle = T,
      callbacks = callback_early_stopping(monitor = "val_loss", 
                                          patience = 15,
                                          restore_best_weights=TRUE,
                                          verbose = verbo)
    )
  
  return(model_localGLMnet)
  
}







# ---------------------------------------

#Encoding Function
preproc = function(
    dt_frame,       # a dataframe
    y,              # target colname required if present in x or for target enc
    num = NULL,     # type of encoding for numericals - normalize, standardize: c("norm","stand") 
    cat = NULL,     # type of encoding for categoricals - one hot, target encoding, entity embedding, joint embedding: c("ohe","targ","entem","joint")
    bypass = NULL,  # string vector of column names to bypass preprocessing
    exclude = NULL, # columns to remove
    verbose = T
){
  
  dt_frame = dt_frame[,which(!(colnames(dt_frame) %in% exclude))]
  
  num_cols = setdiff(dt_frame %>% select_if(is.numeric)  %>% colnames() ,bypass)
  cat_cols = setdiff(colnames(dt_frame),c(num_cols,bypass))
  num_cols = num_cols[num_cols != y]
  
  num_enc_dt = dt_frame %>% select_at(num_cols)
  cat_enc_dt = dt_frame %>% select_at(cat_cols)
  
  
  if(is.null(num)){
    
    if(verbose==T){message("Numerical encoding bypassed")}
    
    num_encoder = function(input){NULL}
    
  }else if(num=="stand"){
    
    num_encoder = function(input){
      
      # ensures statistics are taken from the training dataset 
      stats = apply(num_enc_dt,2,function(x){list(m = mean(x,na.rm = T),s = sqrt(var(x,na.rm = T)))})
      
      toreturn = lapply(num_cols,function(x){data.frame((pull(input,x) - stats[[x]]$m)/stats[[x]]$s) %>% set_names(x)})
      
      return(data.frame(toreturn))
      
    }
    
  }else if(num=="norm"){
    
    num_encoder = function(input){
      
      # ensures statistics are taken from the training dataset 
      stats = apply(num_enc_dt,2,function(x){list(min = min(x,na.rm = T),max = max(x,na.rm = T))})
      toreturn = lapply(num_cols,function(x){data.frame((pull(input,x) - stats[[x]]$min)/(stats[[x]]$max - stats[[x]]$min)) %>% set_names(x)})
      
      return(data.frame(toreturn))
      
    }
    
  }else{
    stop("unrecognized numerical encoding")
  }
  
  
  if(is.null(cat)){
    
    if(verbose==T){message("Categorical encoding bypassed")}
    
    cat_encoder = function(input){NULL}
    
  }else if(cat=="ohe"){
    
    # we have to ensure  levels of categorical vars in test are a subset of  levels in train
    cat_encoder = function(input){
      
      unq_cat = lapply(cat_cols,function(x){sort(unique(cat_enc_dt[,x]))})
      names(unq_cat) = cat_cols
      
      lapply(cat_cols,function(x){
        
        grid = matrix(rep(unq_cat[[x]],length(input[,x])),
                      nrow = length(input[,x]),
                      byrow = T)
        
        dat = matrix(rep(input[,x],length(unq_cat[[x]])),
                     byrow = F,
                     ncol = length(unq_cat[[x]]),
                     dimnames = list(row_names = NULL,col_names = paste0(x,"_",unq_cat[[x]])))
        
        return(data.frame((dat == grid)*1))
        
      }) %>% 
        data.frame() %>% 
        return()
      
    }
    
  }else if(cat=="targ"){
    
    cat_encoder = function(input){
      
      # create lookup
      lkp = apply(cat_enc_dt,
                  2,
                  function(x){
                    
                    data.frame(original = x,
                               target = dt_frame[,y]) %>% 
                      group_by(original) %>% 
                      mutate(enc = mean(target)) %>% 
                      ungroup() %>% 
                      select(-target) %>% 
                      distinct() %>% 
                      arrange(original) %>% 
                      data.frame()
                    
                  })
      
      # replace entries from input according to lkp
      lapply(cat_cols,
             function(x){
               
               plyr::mapvalues(x = input %>% pull(x),
                               from = lkp[[x]]$original,
                               to = lkp[[x]]$enc) %>% 
                 as.numeric() %>% 
                 data.frame() %>% 
                 set_names(x)
               
             }) %>% 
        data.frame() %>% 
        return()
      
    }
    
  }else if(cat=="entem"){
    
  }else if(cat=="joint"){
    
    if(verbose==T){message("Joint embedding is trained...")}
    
    # OHE the categoricals - reference self with categoricals = cat_cols 
    give_ohe = preproc(dt_frame = dt_frame[,cat_cols],y = y,cat = "ohe",num = NULL,verbose = F)
    
    # just the target is going to be present
    ohe_dt = give_ohe(dt_frame)[,-1]
    
    # train the encoder on OHE categoricals
    no_neurons=16
    epoch=200
    batch_size=1000
    learning_rate=0.001
    
    #Network for the autoencoder
    Input = layer_input(shape = c(ncol(ohe_dt)))
    
    Output = Input %>% 
      layer_dense(units=no_neurons, activation='linear', use_bias=FALSE,name="encoder") %>% 
      layer_dense(units=ncol(ohe_dt), activation='softmax', use_bias=TRUE)
    
    model_ae=keras_model(inputs=Input,outputs=Output)
    
    #Optimize the cross entropy
    model_ae %>% 
      compile(optimizer=optimizer_nadam(lr=learning_rate),
              loss="categorical_crossentropy")  %>% 
      fit(ohe_dt,ohe_dt,
          epochs=epoch,
          batch_size=batch_size,
          verbose=0,
          validation_data=list(ohe_dt, ohe_dt),
          callbacks=list(callback_early_stopping(monitor="val_loss", 
                                                 min_delta=0,
                                                 patience=15, 
                                                 verbose=0, 
                                                 mode=c("min"),
                                                 restore_best_weights=TRUE)))
    
    #Recover the representation from the AE
    joint_embedding = keras_model(inputs=model_ae$input, outputs=get_layer(model_ae, "encoder")$output)
    
    cat_encoder = function(input){
      
      joint_embedding %>% 
        predict(give_ohe(input)) %>% 
        return()
      
    }
    
  }else{
    stop("unrecognized categorical encoding")
  }
  
  # combine encoders  
  encoder = function(x){
    
    x = x[,which(!(colnames(x) %in% exclude))]
    
    if(is.null(bypass)){
      bypassed = NULL
    }else if(length(bypass)==1){
      bypassed = x[bypass]
    }else{
      bypassed = x[,bypass]
    }
    
    if(!(y %in% colnames(x))){
      target = NULL
    }else{
      target = x[y]
    }
    
    temp_list = list(target,
                     bypassed,
                     num_encoder(input = x[,num_cols]),
                     cat_encoder(input = x[,cat_cols]))
    
    toreturn = do.call(cbind,
                       temp_list[!unlist(lapply(temp_list,is.null))])
    
    rownames(toreturn) = NULL
    
    return(data.matrix(toreturn))
    
  }
  
  return(encoder)
  
}


# Poisson Deviance 
custom_poisson <- function( y_true, y_pred ) {
  # Mario V. Wüthrich , Michael Merz
  # Statistical Foundations of Actuarial Learning and its Applications
  # Table 4.1 page 87
  
  # 2 * (y_pred - y_true - y_true*log(y_pred/y_true))
  # 2 (μ − y − ylog(μ/y))
  
  K <- backend()
  
  K$mean(2 * (y_pred - y_true - y_true * K$log((y_pred+10^-7)/(y_true+10^-7))))
  
}



# --------------------------------------
encoder = preproc(dt_frame = dt_list$fre_mtpl2_freq[train_rows,],
                  y = "ClaimNb",
                  num = "norm",
                  cat = "ohe",
                  bypass = NULL,
                  exclude = c("IDpol","Exposure"),
                  verbose = T)


CV = 5 #K Folds

CV_vec = sample(1:CV,replace = T,size = nrow(dt_list$fre_mtpl2_freq))
train_rows = which(CV_vec != 1)

train = encoder(dt_list$fre_mtpl2_freq[train_rows,])
test = encoder(dt_list$fre_mtpl2_freq[-train_rows,])

train_LocalGLMnet(
  dt = train[,-1],
  y = train[,1],
  vdt = list(x_val = test[,-1],
             y_val = test[,1]),
  n_dropout = c(F,F,F),
  n_layers = 3,
  n_units = c(20,15,15),
  n_activ = c("tanh","tanh","tanh"),
  lr = 0.005,
  bs = 2^12,
  ep = 1000)
