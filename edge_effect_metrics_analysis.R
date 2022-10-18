
##%######################################################%##
#                                                          #
####            INLA models: Quantifying the            ####
####           drivers of field-edge effects            ####
#                                                          #
##%######################################################%##

# rpr = Relative patch richness
# shdi = Shannon's diversity index
# shei = Shannon's eveness index

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

# Load BP and Beta datasets

bpdf <- as.data.frame(read_csv('./data/data_cleaning/ms_analysis/mcmc_bp_data.csv')) %>% 
  filter(grad >= 0, hole == 'No') %>% 
  mutate(sample_year_f = as.factor(sample_year),
         second_cropping = as.factor(second_cropping),
         py_aes = as.factor(py_aes),
         field_id = field)

bpdf$sample_year_f <- relevel(bpdf$sample_year_f, '2012')
bpdf$py_aes <- relevel(bpdf$py_aes, '0')

betadf <- as.data.frame(read_csv('./data/data_cleaning/ms_analysis/mcmc_beta_data.csv')) %>% 
  filter(hole == 'No') %>% 
  mutate(sample_year_f = as.factor(sample_year),
         second_cropping = as.factor(second_cropping),
         py_aes = as.factor(py_aes),
         field_id = field)

betadf$sample_year_f <- relevel(betadf$sample_year_f, '2012')
betadf$py_aes <- relevel(betadf$py_aes, '0')

###

osr_bp <- droplevels(as.data.frame(bpdf[bpdf$type_of_crop == 'W-OSR',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, bp_dist, bp_dist_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f) %>% 
  drop_na()
osr_beta <- droplevels(as.data.frame(betadf[betadf$type_of_crop == 'W-OSR',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, beta, beta_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f) %>% 
  drop_na()

whe_bp <- droplevels(as.data.frame(bpdf[bpdf$type_of_crop == 'W-Wheat',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, bp_dist, bp_dist_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f, second_cropping) %>% 
  drop_na()
whe_beta <- droplevels(as.data.frame(betadf[betadf$type_of_crop == 'W-Wheat',])) %>% 
  dplyr::select(field, type_of_crop, sample_year, beta, beta_er, py_aes, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi, rainfall_mm, temperature_degc, sample_year_f, second_cropping) %>% 
  drop_na()

##%######################################################%##
#                                                          #
####                  WHEAT: BP models                  ####
#                                                          #
##%######################################################%##

whe_bp$rainfall_mm_p1<- poly(whe_bp$rainfall_mm, 2)[,1]
whe_bp$rainfall_mm_p2 <- poly(whe_bp$rainfall_mm, 2)[,2]
whe_bp$temperature_degc_p1 <- poly(whe_bp$temperature_degc, 2)[,1]
whe_bp$temperature_degc_p2 <- poly(whe_bp$temperature_degc, 2)[,2]

u <- sd(whe_bp$bp_dist)
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(u, 0.01)))

set.seed(3103)
m0 <- inla(bp_dist ~ 1 + 
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

set.seed(3103)
m0b <- inla(bp_dist ~ second_cropping + 
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

set.seed(3103)
m1 <- inla(bp_dist ~ sample_year_f + 
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE, cpo = T),
           family = 'gaussian',
           data = whe_bp)

wht_mod_comp(m0, m0b, m1) 

set.seed(3103)
m1_ar1 <- inla(bp_dist ~
                 f(sample_year, model ='ar1') + 
                 f(field, model = 'iid'), 
               scale = 1/bp_dist_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = whe_bp)

set.seed(3103)
m1_rw1 <- inla(bp_dist ~
                 f(sample_year, model ='rw1', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/bp_dist_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = whe_bp)

wht_mod_comp(m0, m1, m1_ar1, m1_rw1) 

set.seed(3103)
m2a <- inla(bp_dist ~ rainfall_mm_p1 +
             f(sample_year, model ='ar1') +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

set.seed(3103)
m2b <- inla(bp_dist ~ rainfall_mm_p2 +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

set.seed(3103)
m2c <- inla(bp_dist ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.predictor = list(compute = T),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

wht_mod_comp(m1_ar1, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(bp_dist ~ temperature_degc_p1 +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

set.seed(3103)
m3b <- inla(bp_dist ~ temperature_degc_p2 +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

set.seed(3103)
m3c <- inla(bp_dist ~ temperature_degc_p1 +
              temperature_degc_p2 +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

wht_mod_comp(m1_ar1, m3a, m3b, m3c) 

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c) 

set.seed(3103)
m6 <- inla(bp_dist ~ rainfall_mm_p1 +
             rainfall_mm_p2 +
             py_aes +
             f(sample_year, model ='ar1') +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)
  
wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m6) 

set.seed(3103)
m7 <- inla(bp_dist ~ rainfall_mm_p1 +
             rainfall_mm_p2 +
             snh_non_crop +
             f(sample_year, model ='ar1') +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

wht_mod_comp(m2c, m7)

set.seed(3103)
m8 <- inla(bp_dist ~ rainfall_mm_p1 +
             rainfall_mm_p2 +
             snh_all_snh +
             f(sample_year, model ='ar1') +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

wht_mod_comp(m2c, m8) 

set.seed(3103)
m9 <- inla(bp_dist ~ rainfall_mm_p1 +
             rainfall_mm_p2 +
             snh_top_snh +
             f(sample_year, model ='ar1') +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_bp)

wht_mod_comp(m2c, m9)

m10 <- inla(bp_dist ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              rpr +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

wht_mod_comp(m2c, m10) 

m11 <- inla(bp_dist ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              shdi +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

wht_mod_comp(m2c, m11) 

m12 <- inla(bp_dist ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              shei +
              f(sample_year, model ='ar1') +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

wht_mod_comp(m2c, m12) 

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m6, m7, m8, m9, m10, m11, m12)

#
whe_bp %>% 
  select(bp_dist, bp_dist_er, rainfall_mm, temperature_degc, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi) %>% 
  cor(use = 'complete.obs', method = 'pearson')
#

m13 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              shdi +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

m13b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
               py_aes +
               snh_top_snh * shdi +
               sample_year_f +
               f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_bp)

m14 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              shei +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

m14b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
               py_aes +
               snh_top_snh * shei +
               sample_year_f +
               f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_bp)

m15 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              rpr +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_bp)

m15b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p2 +
               py_aes +
               snh_top_snh * rpr +
               sample_year_f +
               f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_bp)

wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m6, m7, m8, m9, m10, m11, m12, m13, m13b, m14, m14b, m15, m15b)

mod_tab <- wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m6, m7, m8, m9, m10, m11, m12, m13, m13b, m14, m14b, m15, m15b)

rm(list=ls(pattern="^m"))

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
                 f(sample_year, model ='ar1') + 
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
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m2b <- inla(beta ~ rainfall_mm_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m2c <- inla(beta ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1_ar1, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(beta ~ temperature_degc_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m3b <- inla(beta ~ temperature_degc_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m3c <- inla(beta ~ temperature_degc_p1 +
              temperature_degc_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1_ar1, m3a, m3b, m3c) 

set.seed(3103)
m4a <- inla(beta ~ rainfall_mm_p1 +
              temperature_degc_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

set.seed(3103)
m4b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b) 

set.seed(3103)
m6 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6) 

set.seed(3103)
m7 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes +
             snh_non_crop +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m7b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes * snh_non_crop +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

wht_mod_comp(m6, m7, m7b) 

set.seed(3103)
m8 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes +
             snh_all_snh +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m8b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes * snh_all_snh +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m8, m8b) 

set.seed(3103)
m9 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes +
             snh_top_snh +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

set.seed(3103)
m9b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
             py_aes * snh_top_snh +
             sample_year_f +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = whe_beta)

wht_mod_comp(m6, m9, m9b) 

m10 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              rpr +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m10b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes * rpr +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m10, m10b) 

m11 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              shdi +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m11b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes * shdi +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

wht_mod_comp(m6, m11, m11b) 

m12 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              shei +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m12b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes * shei +
              sample_year_f +
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

m13 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              shdi +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m13b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh * shdi +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m14 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              shei +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m14b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
               py_aes +
               snh_top_snh * shei +
               sample_year_f +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)

m15 <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
              py_aes +
              snh_top_snh +
              rpr +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = whe_beta)

m15b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 +
               py_aes +
               snh_top_snh * rpr +
               sample_year_f +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = whe_beta)


mod_tab <- wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6, m7, m7b, m8, m8b, m9, m9b, m10, m10b, m11, m11b, m12, m12b, m13, m13b, m14, m14b, m15, m15b)

rm(list=ls(pattern="^m"))

##%######################################################%##
#                                                          #
####                 OILSEED: BP models                 ####
#                                                          #
##%######################################################%##

osr_bp$rainfall_mm_p1<- poly(osr_bp$rainfall_mm, 2)[,1]
osr_bp$rainfall_mm_p2 <- poly(osr_bp$rainfall_mm, 2)[,2]
osr_bp$temperature_degc_p1 <- poly(osr_bp$temperature_degc, 2)[,1]
osr_bp$temperature_degc_p2 <- poly(osr_bp$temperature_degc, 2)[,2]

u <- sd(osr_bp$bp_dist)
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(u, 0.01)))

set.seed(3103)
m0 <- inla(bp_dist ~ 1 + 
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_bp)

set.seed(3103)
m1 <- inla(bp_dist ~ sample_year_f + 
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE, cpo = T),
           family = 'gaussian',
           data = osr_bp)

wht_mod_comp(m0, m1) 

set.seed(3103)
m1_ar1 <- inla(bp_dist ~
                 f(sample_year, model ='ar1') + 
                 f(field, model = 'iid'), 
               scale = 1/bp_dist_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = osr_bp)

set.seed(3103)
m1_ar1p <- inla(bp_dist ~
                 f(sample_year, model ='ar1', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/bp_dist_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = osr_bp)

set.seed(3103)
m1_rw1 <- inla(bp_dist ~
                 f(sample_year, model ='rw1', hyper = hyper.prec) + 
                 f(field, model = 'iid'), 
               scale = 1/bp_dist_er, 
               quantiles = c(0.025, 0.975),
               control.compute=list(dic=TRUE, waic=TRUE),
               family = 'gaussian',
               data = osr_bp)

wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m1_ar1p) 

set.seed(3103)
m2a <- inla(bp_dist ~ rainfall_mm_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m2b <- inla(bp_dist ~ rainfall_mm_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m2c <- inla(bp_dist ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m1_ar1, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(bp_dist ~ temperature_degc_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m3b <- inla(bp_dist ~ temperature_degc_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m3c <- inla(bp_dist ~ temperature_degc_p1 +
              temperature_degc_p2 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m1_ar1, m3a, m3b, m3c) 

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c)

set.seed(3103)
m4a <- inla(bp_dist ~ rainfall_mm_p1 +
              temperature_degc_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m4b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
              sample_year_f +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m4a, m4b)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b) 

set.seed(3103)
m6 <- inla(bp_dist ~ sample_year_f + 
             rainfall_mm_p1 * temperature_degc_p1 +
             py_aes +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_bp)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6) 

set.seed(3103)
m7 <- inla(bp_dist ~ sample_year_f + 
             rainfall_mm_p1 * temperature_degc_p1 +
             snh_non_crop +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_bp)

wht_mod_comp(m4b, m7) 

set.seed(3103)
m8 <- inla(bp_dist ~ sample_year_f + 
             rainfall_mm_p1 * temperature_degc_p1 +
             snh_all_snh +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_bp)

wht_mod_comp(m4b, m8) 

set.seed(3103)
m9 <- inla(bp_dist ~ sample_year_f + 
             rainfall_mm_p1 * temperature_degc_p1 +
             snh_top_snh +
             f(field, model = 'iid'), 
           scale = 1/bp_dist_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_bp)

wht_mod_comp(m4b, m9)

m10 <- inla(bp_dist ~ sample_year_f + 
              rainfall_mm_p1 * temperature_degc_p1 +
              rpr +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m4b, m10) 

m11 <- inla(bp_dist ~ sample_year_f + 
              rainfall_mm_p1 * temperature_degc_p1 +
              shdi +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m4b, m11) 

m12 <- inla(bp_dist ~ sample_year_f + 
              rainfall_mm_p1 * temperature_degc_p1 +
              shei +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

wht_mod_comp(m4b, m12) 

#
osr_bp %>% 
  select(bp_dist, bp_dist_er, rainfall_mm, temperature_degc, snh_non_crop, snh_all_snh, snh_top_snh, rpr, shei, shdi) %>% 
  cor(use = 'complete.obs', method = 'pearson')
#

set.seed(3103)
m13 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
              sample_year_f +
              snh_top_snh +
              shdi +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m13b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
               sample_year_f +
               snh_top_snh * shdi +
             f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_bp)

set.seed(3103)
m14 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
              sample_year_f +
              snh_top_snh +
              shei +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m14b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
               sample_year_f +
               snh_top_snh * shei +
               f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_bp)

set.seed(3103)
m15 <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
              sample_year_f +
              snh_top_snh +
              rpr +
              f(field, model = 'iid'), 
            scale = 1/bp_dist_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_bp)

set.seed(3103)
m15b <- inla(bp_dist ~ rainfall_mm_p1 * temperature_degc_p1 +
               sample_year_f +
               snh_top_snh * rpr +
               f(field, model = 'iid'), 
             scale = 1/bp_dist_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_bp)

wht_mod_comp(m4b,m6, m7, m8, m9, m10, m11, m12, m13, m13b, m14, m14b, m15, m15b) 

mod_tab <- wht_mod_comp(m0, m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m6, m7, m8, m9, m10, m11, m12, m13, m13b, m14, m14b, m15, m15b)


rm(list=ls(pattern="^m"))

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

wht_mod_comp(m0, m1, m1_ar1, m1_ar1p, m1_rw1) 

set.seed(3103)
m2a <- inla(beta ~ rainfall_mm_p1 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m2b <- inla(beta ~ rainfall_mm_p2 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m2c <- inla(beta ~ rainfall_mm_p1 +
              rainfall_mm_p2 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m1_ar1, m2a, m2b, m2c) 

set.seed(3103)
m3a <- inla(beta ~ temperature_degc_p1 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m3b <- inla(beta ~ temperature_degc_p2 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m3c <- inla(beta ~ temperature_degc_p1 +
              temperature_degc_p2 +
              f(sample_year, model = "ar1") +
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
              temperature_degc_p2 +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m4b <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p1 + 
              rainfall_mm_p2 + 
              temperature_degc_p2 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m4c <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 + 
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)


set.seed(3103)
m4d <- inla(beta ~ rainfall_mm_p1 * temperature_degc_p2 + 
              rainfall_mm_p2 + 
              temperature_degc_p1 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)


set.seed(3103)
m4e <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p1 + 
              rainfall_mm_p1 + 
              temperature_degc_p2 +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m4c, m4d, m4e) 

set.seed(3103)
m6 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'), 
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

wht_mod_comp(m1, m1_ar1, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m4c, m4d, m4e, m6) 

set.seed(3103)
m7 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_non_crop +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m7b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
             py_aes * snh_non_crop +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

(x <- wht_mod_comp(m6, m7, m7b)) 

set.seed(3103)
m8 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_all_snh +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m8b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
             py_aes * snh_all_snh +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

wht_mod_comp(m6, m8, m8b) 

set.seed(3103)
m9 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
             rainfall_mm_p1 + 
             temperature_degc_p1 +
             py_aes +
             snh_top_snh +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

set.seed(3103)
m9b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
             py_aes * snh_top_snh +
             f(sample_year, model = "ar1") +
             f(field, model = 'iid'),
           scale = 1/beta_er, 
           quantiles = c(0.025, 0.975),
           control.compute=list(dic=TRUE, waic=TRUE),
           family = 'gaussian',
           data = osr_beta)

wht_mod_comp(m6, m9, m9b)

m10 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              rpr +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

m10b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
              py_aes * rpr +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m6, m10, m10b) 

m11 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              shdi +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

m11b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
              py_aes * shdi +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

wht_mod_comp(m6, m11, m11b)

m12 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              shei +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

m12b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
              py_aes * shei +
              f(sample_year, model = "ar1") +
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
m13 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
                rainfall_mm_p1 + 
                temperature_degc_p1 +
                py_aes +
                snh_top_snh +
                shdi +
                f(sample_year, model = "ar1") +
                f(field, model = 'iid'), 
              scale = 1/beta_er, 
              quantiles = c(0.025, 0.975),
              control.compute=list(dic=TRUE, waic=TRUE),
              family = 'gaussian',
              data = osr_beta)

set.seed(3103)
m13b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * shdi +
               f(sample_year, model = "ar1") +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

set.seed(3103)
m14 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              snh_top_snh +
              shei +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m14b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * shei +
               f(sample_year, model = "ar1") +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

set.seed(3103)
m15 <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
              rainfall_mm_p1 + 
              temperature_degc_p1 +
              py_aes +
              snh_top_snh +
              rpr +
              f(sample_year, model = "ar1") +
              f(field, model = 'iid'), 
            scale = 1/beta_er, 
            quantiles = c(0.025, 0.975),
            control.compute=list(dic=TRUE, waic=TRUE),
            family = 'gaussian',
            data = osr_beta)

set.seed(3103)
m15b <- inla(beta ~ rainfall_mm_p2 * temperature_degc_p2 +
               rainfall_mm_p1 + 
               temperature_degc_p1 +
               py_aes +
               snh_top_snh * rpr +
               f(sample_year, model = "ar1") +
               f(field, model = 'iid'), 
             scale = 1/beta_er, 
             quantiles = c(0.025, 0.975),
             control.compute=list(dic=TRUE, waic=TRUE),
             family = 'gaussian',
             data = osr_beta)

wht_mod_comp(m10, m13, m13b, m14, m14b, m15, m15b)  

whe_beta_best_models <- wht_mod_comp(m0, m1, m1_ar1, m1_ar1p, m1_rw1, m2a, m2b, m2c, m3a, m3b, m3c, m4a, m4b, m4c, m6, m7, m7b, m9, m9b, m10, m10b, m11, m11b, m12, m12b, m13, m13b, m14, m14b, m15, m15b)

rm(list=ls(pattern="^m"))

# End