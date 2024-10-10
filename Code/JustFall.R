# VAST EOF model with seasonal effects
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/

# install.packages('TMB', type = 'source')
# remotes::install_github("james-thorson/VAST")
library(dplyr)
library(tidyr)
library(ggplot2)
library(VAST)
library(tidyverse)
library(FishStatsUtils)
library(sf)

df <- read.csv("Data/Derived/CWM_dataset.csv") %>% 
  st_as_sf(coords = c("decdeg_beglon", "decdeg_beglat"), remove = F) %>%
  st_set_crs("NAD83") %>%
  st_join(ecodata::epu_sf %>% st_make_valid()) %>%
  st_drop_geometry() %>%
  dplyr::select(-survey_area, -Shape_Leng, -Shape_Area) %>%
  pivot_longer(cols = trophic_level:max_obs_length, names_to = "trait", values_to = "value")

trait_list <- c("trophic_level", "offspring_size", "age_maturity", "length_maturity", "l_inf", "k", "max_obs_length") # Drop fecundity due to scale issues


vast_wrapper <- function(n_x = 50){
  ## Seasonal model -----
  working_dir <- here::here("EOF_output/analysis/vast_fall_EOF/")
  
  if(!dir.exists(working_dir)) {
    dir.create(working_dir, recursive  = TRUE)
  }
  #
  ## Attempt to create a log file
  my_log <- file(sprintf("%s/log-%s.txt", working_dir, n_x)) # File name of output log
  
  sink(my_log, append = TRUE, type = "output") # Writing console output to log file
  on.exit(sink(file = NULL), add = TRUE, after = TRUE)
  
  sink(my_log, append = TRUE, type = "message")
  on.exit(sink(file = NULL), add = TRUE, after = TRUE)
  
  
  trait_dat <- df %>%
    rename(year = est_year, 
           lat = decdeg_beglat, 
           lon = decdeg_beglon) %>%
    dplyr::filter(trait %in% trait_list,
                  EPU %in% c("GB", "GOM", "MAB"), 
                  year >=2005, # Do this to increase model convergence speed only...
                  season == "Fall") %>% 
    dplyr::mutate(areaswept_km2 = 1,
                  year_level = factor(year),
                  trait_number = as.numeric(factor(trait)) - 1) %>%
    select(year,
           lat,
           lon,
           areaswept_km2,
           trait_number,
           trait,
           value) %>%
    data.frame()
  
  
  #####
  ## Model settings
  #####
  
  
  ##  Random Fields -----
  
  ## Control the random fields part of the model.
  ## Omega = X is the number of random spatial fields to apply
  ## and Epsilon = X is the number of random spatio-temporal
  ## fields to apply. Omega1 is for the probability of
  ## occurrence, and Omega2 is for the density given occurrence,
  ## similarly for Epsilon.
  
  ## 0 = off
  ## "AR1" = AR1 process
  ## >1 = number of elements in a factor-analysis covariance
  ## "IID" = random effect following an IID distribution
  
  # FieldConfig <- c("Omega1" = "IID", "Epsilon1" = "IID",
  #                  "Omega2" = "IID", "Epsilon2" = "IID")
  
  # Change this according to VAST EOF tutorial...
  # FieldConfig = matrix(c("IID","Identity","IID",2, 0,0,"IID","Identity"), ncol = 2, nrow = 4,
  #                      dimnames = list(c("Omega","Epsilon","Beta","Epsilon_year"),
  #                                      c("Component_1","Component_2")))
  # 
  # FieldConfig = matrix(c(0,"Identity","IID",2, 0,0,"IID","Identity"), ncol = 2, nrow = 4,
  #                      dimnames = list(c("Omega","Epsilon","Beta","Epsilon_year"),
  #                                      c("Component_1","Component_2")))
  ## Autoregressive structure -----
  
  ## Control autoregressive structure for parameters over time
  ## Changing the settings here creates different
  ## autoregressive models for the intercept (Beta) and
  ## spatio-temporal process (Epsilon).
  
  ## 0 = each year is a fixed effect
  ## 1 = random effect
  ## 2 = random walk
  ## 3 = fixed effect that is constant over time
  ## 4 = AR1 process
  
  # RhoConfig <- c("Beta1"    = c(0, 1, 2, 3, 4)[1],
  #                "Beta2"    = c(0, 1, 2, 3, 4)[1],
  #                "Epsilon1" = c(0, 1, 2, 3, 4)[1],
  #                "Epsilon2" = c(0, 1, 2, 3, 4)[1])
  
  RhoConfig=c("Beta1"=1,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0)
  
  ## Correlated overdispersion -----
  
  ## Control correlated overdispersion among categories
  ## for each level of v_i, where eta1 is for encounter
  ## probability, and eta2 is for positive catch rates
  # eta1 = vessel effects on prey encounter rate
  # eta2 = vessel effects on prey weight
  
  ## 0 = off,
  ## "AR1" = AR1 process,
  ## >0 = number of elements in a factor-analysis covariance
  OverdispersionConfig <- c("eta1" = 0,
                            "eta2" = 0)
  
  
  ## Observation model -----
  # Control observation model structure. The first
  # component sets the distribution of the positive
  # distribution component. ?VAST::make_data()
  
  ObsModel <- c("PosDist" = 1,
                "Link"    = 3)
  
  # Make settings
  settings =  make_settings(n_x = n_x,
                            Region = "northwest_atlantic",
                            strata.limits = "EPU",
                            # strata.limits = list('All_areas' = 1:1e5),
                            purpose = "EOF3",
                            n_categories = 2,
                            # ObsModel = ObsModel,
                            RhoConfig = RhoConfig,
                            use_anisotropy = FALSE,
                            #FieldConfig = FieldConfig,
                            bias.correct = FALSE#,
                            # Options = c('treat_nonencounter_as_zero' = TRUE)
  )
  
  settings$epu_to_use <- c("Georges_Bank", "Gulf_of_Maine", "Mid_Atlantic_Bight")
  
  #####
  ## Model fit -- make sure to use new functions
  #####
  
  # Don't use X_contrasts, so that fixed season-year slope isn't confounded with beta term
  fit = fit_model(settings = settings,
                  Lat_i = trait_dat$lat,
                  Lon_i = trait_dat$lon,
                  # t_i = trait_dat$year,
                  t_i = as.numeric(trait_dat$year),
                  c_i = trait_dat$trait_number,
                  b_i = as_units(trait_dat$value, "kg"),
                  a_i = as_units(trait_dat$areaswept_km2, "km^2"),
                  epu_to_use = settings$epu_to_use,
                  newtonsteps = 0,
                  # covariate_data = data.frame(trait_dat, Lat=trait_dat$lat, Lon=trait_dat$lon, Year=as.numeric(trait_dat$year)-1),
                  # X1_formula = ~ season,
                  # X1config_cp = matrix( 2, nrow=length(trait_list), ncol=2 ),
                  # X_contrasts = list(season = contrasts(trait_dat$season, contrasts = FALSE)),
                  getsd = TRUE,
                  Use_REML = FALSE,
                  run_model = TRUE,
                  working_dir = working_dir,
                  optimize_args = list("lower" = -Inf,
                                       "upper" = Inf))
  
  plot(fit, plot_set=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15))
  #fit$tmb_list$Obj$fn( fit$tmb_list$Obj$par )
  
  saveRDS(fit, file = paste0(working_dir, "fit.rds"))
  return(fit)
}



possibly_vast_wrapper <- purrr::quietly(vast_wrapper)

vast_runs <- possibly_vast_wrapper(n_x = 50)


vast_runs$result
  
fit <- readRDS(here::here("EOF_output/analysis/vast_fall_EOF/fit.rds"))

results = plot( fit,
                check_residuals=FALSE,
                working_dir = "EOF_output/analysis/vast_fall_EOF/",
                plot_list = c(11,12,13,14,15),
                category_names = trait_list,
                year_labels = levels(trait_dat$year_season),
                strata_names =  c("Georges Bank", "Gulf of Maine", "Mid-Atlantic Bight"))















