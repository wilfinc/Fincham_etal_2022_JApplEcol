
##%######################################################%##
#                                                          #
####        MiMIn models: intrinsic yield models        ####
#                                                          #
##%######################################################%##

# rpr = Relative patch richness
# shdi = Shannon's diversity index
# shei = Shannon's eveness index

library(dplyr)
library(readr)
library(lme4)
library(MuMIn)

options(na.action = "na.fail") 

# Load the data
data <- as.data.frame(read_csv('./data/data_cleaning/ms_analysis/unscaled_analysis_data_q95_dist_cutoff.csv')) %>% 
  mutate(lcm = as.factor(lcm),
         th_thresh_20 = as.factor(th_thresh_20),
         field_id = as.factor(field_id),
         fieldb_edge = as.factor(fieldb_edge),
         second_cropping = as.factor(second_cropping),
         sample_year = as.factor(sample_year))

# set reference level for analysis to 2012
levels(data$sample_year)
data$sample_year <- relevel(data$sample_year, '2012')

# scale variables
data <- data %>%   
  group_by(type_of_crop) %>% 
  mutate(yield = dry_yield,
         odist = dist,
         dist = scale(dist),
         wetness = scale(wetness),
         shading = scale(shading),
         slope = scale(slope)) %>%  
  ungroup()

##%######################################################%##
#                                                          #
####                   Oilseed Models                   ####
#                                                          #
##%######################################################%##

osr <- data[data$type_of_crop == 'W-OSR',] 

# subset variables
osr <- osr[, c('yield', 'dist', 'slope', 'southness', 'wetness', 'shading', 'th_thresh_20', 'field_id', 'sample_year')]

# Raw polynomials
osr_full_mod_raw <- lmer(yield ~ poly(dist, 2, raw = T) + poly(slope,2, raw = T)*poly(southness,2, raw = T) + poly(wetness,2, raw = T) + poly(shading,2, raw = T) + th_thresh_20 + (1|field_id) + (1|sample_year), data = osr, REML = F, control = lmerControl(optimizer = "bobyqa"))

# Orthogonal polynomials
osr_full_mod_orth <- lmer(yield ~ poly(dist, 2) + poly(slope,2)*poly(southness,2) + poly(wetness,2) + poly(shading,2) + th_thresh_20 + (1|field_id) + (1|sample_year), data = osr, REML = F, control = lmerControl(optimizer = "bobyqa"))

job::job(osr_mumin_out_orth = {osr_mumin_orth <- dredge(osr_full_mod_orth, trace = T)}, import = c('osr', 'osr_full_mod_orth'), title = 'OSR MuMIn model dredge ORTH polys')

results <- osr_mumin_out_orth$osr_mumin_orth # max model
saveRDS(results, './data/data_cleaning/ms_analysis/models/mumin/osr_mumin_orth_polys.rds')


##%######################################################%##
#                                                          #
####                    Wheat Models                    ####
#                                                          #
##%######################################################%##

whe <- data[data$type_of_crop == 'W-Wheat',]

# subset variables
whe <- whe[, c('yield', 'dist', 'second_cropping', 'slope', 'southness', 'wetness', 'shading', 'th_thresh_20', 'field_id', 'sample_year')]
  
# Raw polynomials
whe_full_mod_raw <- lmer(yield ~ poly(dist, 2, raw = T) + second_cropping + poly(slope,2, raw = T) * poly(southness,2, raw = T) + poly(wetness,2, raw = T) + poly(shading,2, raw = T) + th_thresh_20 + (1|field_id) + (1| sample_year), data = whe, REML = F, control = lmerControl(optimizer = "bobyqa"))

# Orthogonal polynomials
whe_full_mod_orth <- lmer(yield ~ poly(dist, 2) + second_cropping + poly(slope,2)*poly(southness,2) + poly(wetness,2) + poly(shading,2) + th_thresh_20 + (1|field_id) + (1| sample_year), data = whe, REML = F, control = lmerControl(optimizer = "bobyqa"))

job::job(wwhe_mumin_out_orth = {wwhe_mumin_orth <- dredge(whe_full_mod_orth, fixed = c('second_cropping'), trace = T)}, import = c('whe', 'whe_full_mod_orth'), title = 'WWHE MuMIn model dredge ORTH polys')

results <- wwhe_mumin_out_orth$wwhe_mumin_orth # max model
saveRDS(results, './data/data_cleaning/ms_analysis/models/mumin/whe_mumin_orth_polys.rds')

# END 