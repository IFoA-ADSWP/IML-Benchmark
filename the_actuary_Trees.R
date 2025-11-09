source("init.R")

CV = 5

set.seed(1)

CV_vec = sample(1:CV,replace = T,size = nrow(dt_list$fre_mtpl2_freq))

models = list()
results = list()
losses = data.frame(CV = paste0("CV_",1:CV),
                    homog = NA,
                    GLM = NA,
                    IBLM = NA,
                    XGB = NA,
                    GAM = NA,
                    EBM = NA)

for (i in 1:CV){
  
  valid_rows = which(CV_vec == i)
  test_i = (i %% 5) + 1
  test_rows = which(CV_vec == test_i)
  train_rows = which(CV_vec != i & CV_vec != test_i)
  
  iter = paste0("CV_",i)
  
  models[[iter]] = list()
  
  results[[iter]] = data.frame(ID = dt_list$fre_mtpl2_freq$IDpol[test_rows],
                               actual = dt_list$fre_mtpl2_freq$ClaimNb[test_rows],
                               GLM = NA,
                               IBLM = NA,
                               XGB = NA,
                               GAM = NA,
                               EBM = NA) %>% 
    mutate(homog = mean(dt_list$fre_mtpl2_freq$ClaimNb[train_rows]))
  
  # homogenous model ------------------------------------------------- 
  
  info_helper(n=paste0(iter," homog"))
  
  losses$homog[i] = poisson_deviance(y_true = results[[iter]]$actual,
                                     y_pred = results[[iter]]$homog)
  
  # GLM ------------------------------------------------- 
  
  info_helper(n=paste0(iter," GLM"))
  
  models[[iter]]$glm_model = glm(formula = ClaimNb~.,
                                 family = poisson,
                                 data = dt_list$fre_mtpl2_freq[train_rows,-c(1,3)])
  
  results[[iter]]$GLM = as.vector(predict(models[[iter]]$glm_model,
                                          dt_list$fre_mtpl2_freq[test_rows,-c(1,3)],type="response"))
  
  losses$GLM[i] = poisson_deviance(y_true = results[[iter]]$actual,
                                   y_pred = results[[iter]]$GLM)
  
  # IBLM  -------------------------------------------
  
  info_helper(n=paste0(iter," IBLM"))
  
  models[[iter]]$IBLM = IBLM::train_iblm_xgb(df_list = list(train = dt_list$fre_mtpl2_freq[train_rows,-c(1,3)],
                                                            validate = dt_list$fre_mtpl2_freq[valid_rows,-c(1,3)]),
                                             family = "poisson",
                                             response_var = "ClaimNb")
  
  results[[iter]]$IBLM = as.vector(predict(models[[iter]]$IBLM,
                                           dt_list$fre_mtpl2_freq[test_rows,-c(1,3)],type="response"))
  
  losses$IBLM[i] = poisson_deviance(y_true = results[[iter]]$actual,
                                    y_pred = results[[iter]]$IBLM)
  
  # XGB ------------------------------------------------- 
  
  info_helper(n=paste0(iter," XGB"))
  
  models[[iter]]$XGB = xgb.train(
    data = xgb.DMatrix(data.matrix(dt_list$fre_mtpl2_freq[train_rows,-c(1,2,3)]),
                       label = dt_list$fre_mtpl2_freq$ClaimNb[train_rows]), 
    watchlist = list(validation = xgb.DMatrix(data.matrix(dt_list$fre_mtpl2_freq[valid_rows,-c(1,2,3)]),
                                              label = dt_list$fre_mtpl2_freq$ClaimNb[valid_rows])),
    params = list(objective = "count:poisson",eval_metric = "poisson-nloglik"),
    nrounds = 1000,
    early_stopping_rounds = 10,
    verbose = 1
  )
  
  results[[iter]]$XGB = as.vector(
    predict(models[[iter]]$XGB,
            xgb.DMatrix(data.matrix(dt_list$fre_mtpl2_freq[test_rows,-c(1,2,3)])),
            type="response")
  )
  
  losses$XGB[i] = poisson_deviance(y_true = results[[iter]]$actual,
                                   y_pred = results[[iter]]$XGB)
  
  # GAM ------------------------------------------------- 
  
  # info_helper(n=paste0(iter," GAM"))
  # models[[iter]]$GAM = NA
  # results[[iter]]$GAM = NA
  # losses$GAM[i] = NA
  
  # EBM ------------------------------------------------- 
  
  # info_helper(n=paste0(iter," EBM"))
  # models[[iter]]$EBM = NA
  # results[[iter]]$EBM = NA
  # losses$EBM[i] = NA
    
  
}

sink(NULL)

# save files
# saveRDS(list(losses = losses,
#              results = results,
#              models = models),file = "The_Actuary_IML.rds")

# temp = readRDS("The_Actuary_IML.rds")
# losses = temp$losses
# results = temp$results

# saveRDS(list(losses = losses,
#              results = results),file = "Results/The_Actuary_IML_wo_models.rds")

# check calibration
bind_rows(results,.id = "id") %>% 
  select(-ID) %>% 
  group_by(id) %>% 
  summarise_all(mean)

# avg deviation from actual %
bind_rows(results,.id = "id") %>% 
  select(-ID) %>% 
  group_by(id) %>% 
  summarise_all(mean) %>% 
  mutate_if(is.numeric,~if_else(.==actual,.,./actual - 1)) %>% 
  select(-actual) %>% 
  mutate_if(is.numeric,scales::percent,0.1)

analysis = bind_rows(results,.id = "id")  %>% 
  # select(id,actual,glm,XGB, homog, train_GLM_w_XGB, GLM_XGB,IBLM) %>% 
  pivot_longer(cols = GLM:EBM) %>% 
  filter(!is.na(value)) %>% 
  mutate(actual = actual,
         value = value,
         poiss = Vectorize(poisson_deviance)(y_true = actual,
                                             y_pred = value)) 

# ovarall and per fold results
poiss_per_CV = rbind(losses,
                     losses %>%
                       pivot_longer(cols = !CV) %>%
                       group_by(name) %>%
                       summarise(mean_poiss = mean(value)) %>%
                       arrange(mean_poiss) %>%
                       pivot_wider(values_from = mean_poiss,names_from = name) %>%
                       mutate(CV = "mean_poiss"))

# pinball
poiss_per_CV %>% 
  mutate_if(is.numeric,~ if_else(. == homog, .,1 -  ./homog)) %>% 
  select(-homog) %>% 
  mutate_if(is.numeric,scales::percent,0.1)


losses %>% 
  mutate_if(is.numeric,.funs = function(x)(x*c(data.frame(k=CV_vec) %>% 
                                                 count(k) %>% pull(n)))) %>% 
  janitor::adorn_totals()

analysis %>%
  filter(name!="homog") %>%
  rename(model=name) %>% 
  ggplot(aes(x = poiss,fill=model,color=model,linetype=model))+
  geom_density(alpha=0.3,size=1)+
  ggplot2::scale_fill_manual(values = c("blue","yellow","green","red","white"))+
  xlim(0,0.75)+
  # facet_wrap(~name)+
  # ggdark::dark_theme_classic()+
  # theme(panel.grid.minor = element_line(colour="darkgrey", size=0.01,linetype = 3))+
  ggtitle("Poisson deviance per observation, per model")+
  xlab("Poisson deviance")

analysis %>%
  filter(name!="homog") %>%
  rename(model=name) %>% 
  ggplot(aes(x = poiss))+
  facet_wrap(~model,ncol = 1)+
  geom_density(alpha=0.3,size=0.5)+
  ggplot2::scale_fill_manual(values = c("blue","yellow","green","red","white"))+
  xlim(0,0.75)+
  # facet_wrap(~name)+
  # ggdark::dark_theme_classic()+
  # theme(panel.grid.minor = element_line(colour="darkgrey", size=0.01,linetype = 3))+
  ggtitle("Poisson deviance per observation, per model")+
  xlab("Poisson deviance")

# lift chart
multiple_lift(y_true = bind_rows(results,.id = "id") %>% pull(actual),
              y_pred_df = bind_rows(results,.id = "id") %>% select(GLM,
                                                                   XGB,
                                                                   homog,
                                                                   IBLM))+
  ggtitle("Combined lift chart")+
  xlab("Tiles")+
  ylab("Implied frequency")


analysis %>% 
  slice_sample(prop = 0.25) %>% 
  ggplot(aes(x = actual,y=poiss))+geom_point()+facet_wrap(~name)

