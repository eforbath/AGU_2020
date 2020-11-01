#############################################################
#############################################################
######### AGU 2020 - CYN NDVI and Plant Composition #########
#############################################################
#############################################################


getwd()
setwd("/Users/elenaforbath/Downloads/loranty_lab/CYN_plant_comp_data")

## install packages 
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tiff")
install.packages("rtiff")
install.packages("raster") 
install.packages("sp")
install.packages("rgdal")
install.packages("tidyr")
install.packages("stringr")
install.packages("lme4")
install.packages("cowplot")
install.packages("lsmeans")
install.packages("DHARMa")

library(tiff)
library(rtiff)
library(raster)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(stringr)
library(lme4)
library(cowplot)
library(lme4)
library(lsmeans)
library(DHARMa)

##### raster of site #####
FL016 <- raster("CYN_TR1_FL016M/FL016.tif")  
FL016

FL020 <- raster("CYN_TR1_FL020M/FL020.tif")
FL020

##### percent cover ######
percent_cover <- na.omit(read.csv("percent_cover.csv"))
names(percent_cover)[names(percent_cover) == "Plot.ID"] <- "plot"
names(percent_cover)[names(percent_cover) == "Treatment"] <- "treatment"
names(percent_cover)[names(percent_cover) == "Functional.group"] <- "functional_group"
names(percent_cover)[names(percent_cover) == "percent.cover"] <- "percent_cover"




###### point intercept #####
pt_int<- na.omit(read.csv("pt_intercept.csv"))




##### correlation plots #####
install.packages("ggcorrplot")
install.packages("GGally")
install.packages("tidyr")

library(tidyr)
library(ggcorrplot)
library(GGally)




veg = plant_comp %>% 
  dplyr::select(plot:percent_composition) %>% 
  pivot_wider(names_from = functional_groups, 
              values_from = percent_composition,
              values_fill = 0)

veg

names(veg)


veg_var = veg %>%
  dplyr::select(plot, FL016_ndvi, CON:FORB)
veg_var


corr <- round(cor(veg_var), 1)
corr

ggcorrplot(corr)


##### extracting all values from each plots (distribution of values) ######
detach("package:tidyr", unload = TRUE) # extract wont ork with tidyr


library(sf)

plots_feature <- st_read("plots_feature2/plots_feature2.shp")


st_geometry_type(plots_feature)

plots_feature

ggplot() + 
  geom_sf(data = plots_feature, size = 3, color = "black", fill = "cyan1") + 
  ggtitle("Plot Boundaries") + 
  coord_sf()


## extrat ndvi by feature class

pre <- extract(FL016, plots_feature,
                   small = TRUE,
                   df = TRUE, 
                   factor = TRUE)
names(ndvi_FL016)[names(ndvi_FL016) == "FL016"] <- "FL016_ndvi"
pre <- as.data.frame(pre)


post <- extract(FL020, plots_feature,
                    small = TRUE,
                    df = TRUE, 
                    factor = TRUE)
names(ndvi_FL020)[names(ndvi_FL020) == "FL020"] <- "FL020_ndvi"
post <- as.data.frame(post)


plots_feature <- as.data.frame(plots_feature)
plots_feature$ID <- c(1:23)

all_new <- merge(pre, plots_feature, by = "ID")


boxplot(FL016 ~ plot_num, data = all_new)
