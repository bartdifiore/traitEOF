library(tinyVAST)
library(sf)
library(rnaturalearth)
library(tidyverse)

trait_list <- c("trophic_level", "offspring_size", "age_maturity", "length_maturity", "l_inf", "k", "max_obs_length") # Drop fecundity due to scale issues

# trait_list <- c("length_maturity", "k", "l_inf", "offspring_size") # Drop fecundity due to scale issues

# project data
trait_dat <- read.csv("Data/CWM_dataset.csv") %>% 
  st_as_sf(coords = c("decdeg_beglon", "decdeg_beglat"), remove = F) %>%
  st_set_crs("NAD83") %>%
  st_join(ecodata::epu_sf %>% st_make_valid()) %>%
  st_drop_geometry() %>%
  dplyr::select(-survey_area, -Shape_Leng, -Shape_Area) %>%
  pivot_longer(cols = trophic_level:max_obs_length, names_to = "trait", values_to = "value") %>%
  rename(year = est_year, 
         lat = decdeg_beglat, 
         lon = decdeg_beglon) %>%
  dplyr::filter(trait %in% trait_list,
                EPU %in% c("GB", "GOM", "MAB"), 
                #year >=2005), # Do this to increase model convergence speed only...
  ) %>%
  dplyr::mutate(season = str_to_lower(season), 
                areaswept_km2 = 1,
                year_level = factor(year),
                season = factor(season, levels = c("spring", "fall")),
                year_labels = factor(paste(year, season, sep = "_")),
                year_season = factor(year_labels, levels = paste(rep(levels(year_level),
                                                                     each = nlevels(season)),
                                                                 levels(season),
                                                                 sep="_")),
                trait_number = as.numeric(factor(trait)) - 1) %>%
  filter(season == "fall") %>%
  select(year,
         lat,
         lon,
         trait,
         value) %>%
                data.frame()
  
n_eof = 2
dsem = make_eof_ram( variables = trait_list,
                     times = sort(unique(trait_dat[,'year'])),
                     n_eof = 2,
                     standard_deviations = 0 )

mesh = fmesher::fm_mesh_2d( trait_dat[,c('lon','lat')], cutoff=0.25)
plot(mesh)

# fit model
out = tinyVAST( dsem = dsem,
                sem = "",
                data = trait_dat,
                formula = value ~ 1,
                family = Gamma(link = "log"),
                spatial_graph = mesh,
                space_column = c("lon","lat"),
                variable_column = "trait",
                time_column = "year",
                # distribution_column = "dist",
                times = c(paste0("EOF_",seq_len(n_eof)), sort(unique(trait_dat[,'year']))),
                control = tinyVASTcontrol( profile="alpha_j",
                                         gmrf_parameterization="projection") )

write_rds(out, "Data/EOF_output/all_years.rds")


df_sf <- read.csv("Data/CWM_dataset.csv") %>% 
  st_as_sf(coords = c("decdeg_beglon", "decdeg_beglat"), remove = F) %>%
  st_set_crs("NAD83") %>%
  st_join(ecodata::epu_sf %>% st_make_valid()) 

# Country shapefiles for plotting
sf_maps = ne_countries( return="sf", scale="medium", continent=c("north america") )
sf_maps = st_transform( sf_maps, crs=st_crs(df_sf) )
sf_maps = st_union( sf_maps )


# Shapefile for water
sf_water = st_difference( st_as_sfc(st_bbox(df_sf)), sf_maps )
sf_land = st_difference(st_as_sfc(st_bbox(df_sf)), sf_water)
plot(sf_water)
plot(sf_land)


# Create extrapolation grid
cellsize = 0.1
sf_grid = st_make_grid( df_sf, cellsize=cellsize )
# Restrict to water
grid_i = st_intersects( sf_water, sf_grid )
sf_grid = sf_grid[ unique(unlist(grid_i)) ]

newdata = data.frame( lon = rep(st_coordinates(st_centroid(sf_grid))[,1], times = length(trait_list)), lat = rep(st_coordinates(st_centroid(sf_grid))[,2]), trait = rep(trait_list, each = length(st_coordinates(st_centroid(sf_grid))[,1]))) %>% filter(trait %in% c("length_maturity"))

# newdata <- expand.grid(newdata, trait = sort(unique(trait_dat[,'trait'])))


# Extract loadings
L_tf = matrix( 0, nrow=length(unique(trait_dat$year)), ncol=2,
               dimnames=list(unique(trait_dat$year),c("EOF_1","EOF_2")) )
L_tf[lower.tri(L_tf,diag=TRUE)] = out$opt$par[names(out$opt$par)=="beta_z"]


# Extract factor-responses
EOF1_g = predict( out, cbind(newdata, year="EOF_1"), what="pepsilon_g" )
EOF2_g = predict( out, cbind(newdata, year="EOF_2"), what="pepsilon_g" )
omega_g = predict( out, cbind(newdata, year="EOF_2"), what="pomega_g" )


# Rotate responses and loadings
rotated_results = rotate_pca( L_tf=L_tf, x_sf=cbind(EOF1_g,EOF2_g), order="decreasing" )
#> Warning in sqrt(Eigen$values): NaNs produced
EOF1_g = rotated_results$x_sf[,1]
EOF2_g = rotated_results$x_sf[,2]
L_tf = rotated_results$L_tf


# Plot on map
sf_plot = st_sf( sf_grid, "EOF1_g"=EOF1_g, "EOF2_g"=EOF2_g, "omega_g"=omega_g)

p1 <- ggplot()+
  geom_sf(data = sf_plot, aes(fill = EOF1_g))+
  geom_sf(data = sf_land)+
  scale_fill_viridis_c()+
  theme_bw()

p2 <- ggplot()+
  geom_sf(data = sf_plot, aes(fill = EOF2_g))+
  geom_sf(data = sf_land)+
  scale_fill_viridis_c()+
  theme_bw()

p3 <- ggplot()+
  geom_sf(data = sf_plot, aes(fill = omega_g))+
  geom_sf(data = sf_land)+
  scale_fill_viridis_c()+
  theme_bw()

p4 <- data.frame(year = unique(trait_dat$year), EOF1 = L_tf[,1], EOF2 = L_tf[,2]) %>%
  pivot_longer(EOF1:EOF2)%>%
  ggplot()+
  geom_line(aes(color = name, x = year, y = value))+
  scale_color_manual(values = viridisLite::viridis(n_eof))+
  theme_bw()

cowplot::plot_grid(p1, p2, p3, p4)
ggsave("Figures/EOF_tinyVAST_result1.png", bg = "white")









# 
# par(mfrow=c(2,2), oma=c(2,2,0,0) )
# plot( sf_plot[,'EOF1_g'], reset=FALSE, key.pos=NULL, border=NA )
# plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" )
# plot( sf_plot[,'EOF2_g'], reset=FALSE, key.pos=NULL, border=NA )
# plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" )
# plot( sf_plot[,'omega_g'], reset=FALSE, key.pos=NULL, border=NA )
# plot( st_geometry(sf_maps), add=TRUE, border=NA, col="grey" )
# matplot( y=L_tf, x=unique(trait_dat$year), type="l",
#          col=viridisLite::viridis(n_eof), lwd=2, lty="solid" )
# legend( "top", ncol=n_eof, legend=paste0("EOF",1:n_eof),
#         fill=viridisLite::viridis(n_eof) )
# 
# 









