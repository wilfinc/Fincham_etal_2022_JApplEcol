
############################################################
#                                                          #
#            INLA Models: Sensitivity analyses             #
#                                                          #
############################################################

# rpr = Relative patch richness
# shdi = Shannon's diversity index
# shei = Shannon's eveness index

# Re-running the previous INLA analyses but do not exclude the negative Betas
# BP models will be unchanged

# Functions:

wht_mod_comp <- function(..., metric = 'waic'){ # compares INLA models using WAIC values
  
  model_list <- list(...)
  
  if(metric == 'waic'){
    out <- data.table::rbindlist(
      lapply(model_list, function(x){
        data.frame(model = as.character(x$.args$formula[3])
                   , waic = x$waic$waic, best_model = '.')
      })
    )
    
    out$best_model[which.min(out$waic)] <- '<---'
    
    out <- out %>% arrange(waic)
    
  }
  
  if(metric == 'dic'){
    out <- data.table::rbindlist(
      lapply(model_list, function(x){
        data.frame(model = as.character(x$.args$formula[3])
                   , dic = x$dic$dic, best_model = '.')
      })
    )
    
    out$best_model[which.min(out$dic)] <- '<---'
    
    out <- out %>% arrange(dic)
    
  }
  
  if((sum(sapply(model_list, function(z){length(z$summary.random)})) == length(model_list)) == F | all(c(names(model_list$summary.random), names(model_list$summary.random), names(model_list$summary.random)) == names(model_list$summary.random), na.rm = F) == F){warning('Check that random effects are consistant across models.')}
  
  
  return(out)
}


library(INLA)
library(tidyr)
library(dplyr)
library(readr)

options(scipen = 999)

# Load Beta dataset

betadf <- as.data.frame(read_csv('./data/data_cleaning/ms_analysis/mcmc_beta_data_all.csv')) %>% 
  filter(hole == 'No') %>% 
  mutate(sample_year_f = as.factor(sample_year),
         second_cropping = as.factor(second_cropping),
         py_aes = as.factor(py_aes),
         field_id = field)

betadf$sample_year_f <- relevel(betadf$sample_year_f, '2012')
betadf$py_aes <- relevel(betadf$py_aes, '0')

##

osr_beta <- droplevels(as.data.frame(betadf[betadf$type_of_crop == 'W-OSR',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, beta, beta_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f) %>% 
  drop_na()

whe_beta <- droplevels(as.data.frame(betadf[betadf$type_of_crop == 'W-Wheat',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, beta, beta_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f, second_cropping) %>% 
  drop_na()


##%######################################################%##
#                                                          #
####                 WHEAT: Beta models                 ####
#                                                          #
##%######################################################%##

whe_beta$rainfall_mm_p1<- poly(whe_beta$rainfall_mm, 2)[,1]
whe_beta$rainfall_mm_p2 <- poly(whe_beta$rainfall_mm, 2)[,2]
whe_beta$temperature_degc_p1 <- poly(whe_beta$temperature_degc, 2)[,1]
whe_beta$temperature_degc_p2 <- poly(whe_beta$temperature_degc, 2)[,2]

u <- sd(whe_beta$beta)
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(u, 0.01)))

set.seed(3103)
m0 <- inla(beta ~ 1 + 
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m0b <- inla(beta ~ second_cropping + 
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)


set.seed(3103)
m1 <- inla(beta ~ sample_year_f + 
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE, cpo = T),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m1b <- inla(beta ~ sample_year_f + 
              second_cropping +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE, cpo = T),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m0, m0b, m1, m1b) 

set.seed(3103)
m1_ar1 <- inla(beta ~
                 f(sample_year, model ='ar1',  hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/beta_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = whe_beta)

set.seed(3103)
m1_ar1p <- inla(beta ~
                  f(sample_year, model ='ar1', hyper = hyper.prec) + 
                  f(field, model = 'iid'), 
                scale = 1/beta_er, 
                quantiles = c(0.025, 0.975),
                control.compute=list(dic=TRUE, waic=TRUE),
                family = 'gaussian',
                data = whe_beta)

set.seed(3103)
m1_rw1 <- inla(beta ~
                 f(sample_year, model ='rw1', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/beta_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = whe_beta)

set.seed(3103)
m1_rw2 <- inla(beta ~
                 f(sample_year, model ='rw2', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/beta_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = whe_beta)

wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m1_rw2, m1_ar1p) 

set.seed(3103)
m2a <- inla(beta ~ rainfall_mm_p1 +
              f(sample_year, model ='ar1', hyper = hyper.prec) + 
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m2b <- inla(beta ~ rainfall_mm_p2 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m2c <- inla(beta ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1_ar1p, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(beta ~ temperature_degc_p1 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m3b <- inla(beta ~ temperature_degc_p2 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m3c <- inla(beta ~ temperature_degc_p1 +
              temperature_degc_p1 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1_ar1p, m3a, m3b, m3c) 

set.seed(3103)
m4a <- inla(beta ~ rainfall_mm_p1 +
              temperature_degc_p2 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m4b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1, m1_ar1p, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b) 

set.seed(3103)
m6 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
             py_aes +
             f(sample_year, model ='ar1', hyper = hyper.prec) +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

wht_mod_comp(m1, m1_ar1p, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6) 

wht_mod_comp(m1_ar1p, m4a, m4b, m6) 

set.seed(3103)
m7 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
             py_aes +
             snh_non_crop +
             f(sample_year, model ='ar1', hyper = hyper.prec) +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m7b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes * snh_non_crop +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m7, m7b)

set.seed(3103)
m8 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
             py_aes +
             snh_all_snh +
             f(sample_year, model ='ar1', hyper = hyper.prec) +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m8b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes * snh_all_snh +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m8, m8b) 

set.seed(3103)
m9 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
             py_aes +
             snh_top_snh +
             f(sample_year, model ='ar1', hyper = hyper.prec) +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m9b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes * snh_top_snh +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m9, m9b) 

set.seed(3103)
m10 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes +
              rpr +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m10b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * rpr +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

wht_mod_comp(m6, m10, m10b) 

set.seed(3103)
m11 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes +
              shdi +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m11b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * shdi +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

wht_mod_comp(m6, m11, m11b) 

set.seed(3103)
m12 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes +
              shei +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m12b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * shei +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

wht_mod_comp(m6, m12, m12b) 

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6, m7, m7b, m8, m8b, m9, m9b, m10, m10b, m11, m11b, m12, m12b)

#
whe_beta %>% 
  select(beta, beta_er, rainfall_mm, temperature_degc, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi) %>% 
  cor(use = 'complete.obs', method = 'pearson')
#

set.seed(3103)
m13 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes * snh_top_snh +
              shdi +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m13b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * snh_top_snh +
               snh_top_snh * shdi +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             control.predictor = list(compute = TRUE, link = 1),
             family = 'gaussian',
             data = whe_beta)

set.seed(3103)
m14 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes * snh_top_snh +
              shei +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m14b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * snh_top_snh +
               snh_top_snh * shei +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

set.seed(3103)
m15 <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
              py_aes * snh_top_snh +
              rpr +
              f(sample_year, model ='ar1', hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m15b <- inla(beta ~ rainfall_mm_p1 + temperature_degc_p2 +
               py_aes * snh_top_snh +
               snh_top_snh * rpr +
               f(sample_year, model ='ar1', hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

wht_mod_comp(m11, m13, m13b, m14, m14b, m15, m15b)

mod_tab <- wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6, m7, m7b, m8, m8b, m9, m9b, m10, m10b, m11, m11b, m12, m12b, m13, m13b, m14, m14b, m15, m15b)

##%######################################################%##
#                                                          #
####                OILSEED: Beta models                ####
#                                                          #
##%######################################################%##

osr_beta$rainfall_mm_p1<- poly(osr_beta$rainfall_mm, 2)[,1]
osr_beta$rainfall_mm_p2 <- poly(osr_beta$rainfall_mm, 2)[,2]
osr_beta$temperature_degc_p1 <- poly(osr_beta$temperature_degc, 2)[,1]
osr_beta$temperature_degc_p2 <- poly(osr_beta$temperature_degc, 2)[,2]

u <- sd(osr_beta$beta)
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(u, 0.01)))

set.seed(3103)
m0 <- inla(beta ~ 1 + 
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m1 <- inla(beta ~ sample_year_f + 
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE, cpo = T),
           family = 'gaussian',
           data = osr_beta)

wht_mod_comp(m0, m1)

set.seed(3103)
m1_ar1 <- inla(beta ~
                 f(sample_year, model ='ar1') + 
                 f(field, model = 'iid'), 
               scale = 1/beta_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = osr_beta)

set.seed(3103)
m1_ar1p <- inla(beta ~
                  f(sample_year, model ='ar1', hyper = hyper.prec) + 
                  f(field, model = 'iid'), 
                scale = 1/beta_er, 
                quantiles = c(0.025, 0.975),
                control.compute=list(dic=TRUE, waic=TRUE),
                family = 'gaussian',
                data = osr_beta)

set.seed(3103)
m1_rw1 <- inla(beta ~
                 f(sample_year, model ='rw1', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/beta_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = osr_beta)

wht_mod_comp(m0, m1, m1_ar1, m1_ar1p, m1_rw1) #

set.seed(3103)
m2a <- inla(beta ~ rainfall_mm_p1 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m2b <- inla(beta ~ rainfall_mm_p2 +
              f(sample_year, model = "ar1",  hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m2c <- inla(beta ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m1_ar1p, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(beta ~ temperature_degc_p1 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m3b <- inla(beta ~ temperature_degc_p2 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m3c <- inla(beta ~ temperature_degc_p1 +
              temperature_degc_p2 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m1_ar1, m3a, m3b, m3c)

set.seed(3103)
m4a <- inla(beta ~ rainfall_mm_p1 +
              rainfall_mm_p2 + 
              temperature_degc_p1 +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m4b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p1 + 
              rainfall_mm_p2 + 
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m4c <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p1 + 
              rainfall_mm_p1 + 
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m4c) 

set.seed(3103)
m6 <- inla(beta ~ rainfall_mm_p2 + 
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             f(sample_year, model = "ar1", hyper = hyper.prec) +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m6) 

set.seed(3103)
m7 <- inla(beta ~ rainfall_mm_p2 + 
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_non_crop +
             f(sample_year, model = "ar1", hyper = hyper.prec) +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m7b <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes * snh_non_crop +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'),
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

(x <- wht_mod_comp(m6, m7, m7b)) 

set.seed(3103)
m8 <- inla(beta ~ rainfall_mm_p2 + 
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_all_snh +
             f(sample_year, model = "ar1", hyper = hyper.prec) +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m8b <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes * snh_all_snh +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'),
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m6, m8, m8b) 

set.seed(3103)
m9 <- inla(beta ~ rainfall_mm_p2 + 
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_top_snh +
             f(sample_year, model = "ar1", hyper = hyper.prec) +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

# set.seed(3103)
# m9b <- inla(beta ~ rainfall_mm_p2 + 
#               rainfall_mm_p1 + 
#               temperature_degc_p1 +
#               py_aes * snh_top_snh +
#               f(sample_year, model = "ar1", hyper = hyper.prec) +
#               f(field, model = 'iid'),
#             scale = 1/beta_er, 
#             quantiles = c(0.025, 0.975),
#             control.compute=list(dic=TRUE, waic=TRUE),
#             family = 'gaussian',
#             data = osr_beta)

set.seed(3103)
m9b <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes * snh_top_snh +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'),
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            control.inla=list(h=0.001),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m6, m9, m9b)

set.seed(3103)
m10 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              rpr +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m10b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes * rpr +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

wht_mod_comp(m6, m10, m10b) 

set.seed(3103)
m11 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              shdi +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m11b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes * shdi +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

wht_mod_comp(m6, m11, m11b) 

set.seed(3103)
m12 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              shei +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m12b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes * shei +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

wht_mod_comp(m6, m12, m12b) 

#
osr_beta %>% 
  select(beta, beta_er, rainfall_mm, temperature_degc, snh_top_snh, snh_non_crop, snh_all_snh, rpr, shei, shdi) %>% 
  cor(use = 'complete.obs', method = 'pearson')
#

set.seed(3103)
m13 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              snh_top_snh +
              shdi +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m13b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * shdi +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

set.seed(3103)
m14 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              snh_top_snh +
              shei +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

# set.seed(3103)
# m14b <- inla(beta ~ rainfall_mm_p2 + 
#                rainfall_mm_p1 + 
#                temperature_degc_p1 +
#                py_aes +
#                snh_top_snh * shei +
#                f(sample_year, model = "ar1", hyper = hyper.prec) +
#                f(field, model = 'iid'), 
#              scale = 1/beta_er, 
#              quantiles = c(0.025, 0.975),
#              control.compute=list(dic=TRUE, waic=TRUE),
#              family = 'gaussian',
#              data = osr_beta) 

set.seed(3103)
m14b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * shei +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             control.inla = list(h = 0.001),
             family = 'gaussian',
             data = osr_beta)

set.seed(3103)
m15 <- inla(beta ~ rainfall_mm_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              snh_top_snh +
              rpr +
              f(sample_year, model = "ar1", hyper = hyper.prec) +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m15b <- inla(beta ~ rainfall_mm_p2 + 
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * rpr +
               f(sample_year, model = "ar1", hyper = hyper.prec) +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

wht_mod_comp(m9, m13, m13b, m14, m14b, m15, m15b) 

osr_beta_best_models <- wht_mod_comp(m0, m1, m1_ar1, m1_ar1p, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m4c, m6, m7, m7b, m9, m9b, m10, m10b, m11, m11b, m12, m12b, m13, m13b, m14, m14b, m14bb, m15, m15b)

# End