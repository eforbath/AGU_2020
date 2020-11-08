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
names(pt_int)[names(pt_int) == "plot"] <- "plot_num"



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





###### extract MEAN ndvi by feature class ######
pre_mean<- extract(FL016, plots_feature,
                      fun = mean,
                      df = TRUE)


post_mean<- extract(FL020, plots_feature,
                   fun = mean,
                   df = TRUE)


plots_feature <- as.data.frame(plots_feature)
plots_feature$ID <- c(1:23)
plots = plots_feature %>% select(plot_num, ID)


pre_mean <- merge (pre_mean, plots, by = "ID")
post_mean <- merge (post_mean, plots, by = "ID")


treatments <- read.csv("treatments.csv")
treatments 
treatments$plot = gsub("P", "", treatments$plot)
names(treatments)[names(treatments) == "plot"] <- "plot_num"


pre_mean <- merge(pre_mean, treatments, by = "plot_num")
post_mean <- merge(post_mean, treatments, by = "plot_num")

all_mean <- merge(pre_mean, post_mean, by = "plot_num")
all_pt_int <-merge(all_mean, pt_int, by = "plot_num")

all_pt_int$ndvi_diff <- all_pt_int$FL016 - all_pt_int$FL020



## boxplot of before and after NDVI
plot(FL016 ~  treatment.x, data = all_mean)
plot(FL020 ~ treatment.x, data = all_mean)


lm(FL016 ~  FL020 + treatment.x, data = all_mean)



## bargraph of before and after clipping
barplot(cbind(FL016, FL020) ~ plot_num, 
        all_mean,
        beside = TRUE,
        las=2)

CT <- subset(all_mean, treatment.x == "CT")
SH <- subset(all_mean, treatment.x == "SH")
GR <- subset(all_mean, treatment.x == "GR")
GS <- subset(all_mean, treatment.x == "GS")

barplot(cbind(FL016, FL020) ~ plot_num, 
        CT,
        beside = TRUE,
        las=2,
        col= c("black", "grey"),
        ylim = c(0,1),
        xlab = "plots",
        ylab = "NDVI", 
        main = "Control plots")
legend("topright", legend = c("pre-clipping", "post-clipping"), fill = c("black", "grey"))

barplot(cbind(FL016, FL020) ~ plot_num, 
        SH,
        beside = TRUE,
        las=2, 
        col= c("black", "grey"),
        ylim = c(0,1), 
        xlab = "plots",
        ylab = "NDVI", 
        main = "Shrub removal plots")
legend("topright", legend = c("pre-clipping", "post-clipping"), fill = c("black", "grey"))

barplot(cbind(FL016, FL020) ~ plot_num, 
        GR,
        beside = TRUE,
        las=2,
        col= c("black", "grey"),
        ylim = c(0,1), 
        xlab = "plots",
        ylab = "NDVI", 
        main = "Grass removal plots")
legend("topright", legend = c("pre-clipping", "post-clipping"), fill = c("black", "grey"))

barplot(cbind(FL016, FL020) ~ plot_num, 
        GS,
        beside = TRUE,
        las=2,
        col= c("black", "grey"),
        ylim = c(0,1), 
        xlab = "plots",
        ylab = "NDVI", 
        main = "Grass and shrub plots") 
legend("topright", legend = c("pre-clipping", "post-clipping"), fill = c("black", "grey"))


####all_mean %>% 
####  mutate(ndvi_diff = FL020 - FL016)

lm_mean <- lm(FL020 ~ treatment, data = post_mean)
summary(lm_mean)


##### correlation plots #####
install.packages("ggcorrplot")
install.packages("GGally")
install.packages("tidyr")

library(tidyr)
library(ggcorrplot)
library(GGally)




veg = pt_int %>% 
  dplyr::select(plot_num:percent_composition) %>% 
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


##### extract ALL ndvi by feature class ####

pre <- extract(FL016, plots_feature,
                   small = TRUE,
                   df = TRUE, 
                   factor = TRUE)
pre <- as.data.frame(pre)


post <- extract(FL020, plots_feature,
                    small = TRUE,
                    df = TRUE, 
                    factor = TRUE)
post <- as.data.frame(post)



FL016 <- merge(pre, plots, by = "ID")
FL020 <- merge(post, plots, by = "ID")

hist(subset(FL016$FL016, FL016$plot_num == "102"))
hist(subset(FL016$FL016, FL016$plot_num == "103"))



write.csv(FL016, "FL016.csv") ## have to manually create dataframe
write.csv(FL020, "FL020.csv")

FL016_new <- read.csv("FL016_new.csv")
FL020_new <- read.csv("FL020_new.csv")



boxplot(FL016 ~ plot_num, data = all_new) 



## for loop to create histogram of NDVI values for each plot 

install.packages("reshape") 
library(reshape) 

pre2 <-unstack(pre, form = formula(FL016 ~ ID), drop = TRUE) 


pre3 <- as.data.frame(pre2,row.names = FALSE) ## Error in (function (..., row.names = NULL, check.rows = FALSE, check.names = TRUE,  : 
                            ##                    arguments imply differing number of rows
pre3 <- table(pre2, useNA= "ifany")



for (i in 1:length(pre2)) { # for every column in the "new" data frame
  x <- pre2[ ,i] # identifying columns (?)
  # Plot histogram of x
  jpeg(file = paste("dist", names((pre2)[i]), ".jpeg", sep = ""))
  hist(x,
       main = paste("NDVI distribution for P", names((pre2)[i])), #paste name of column to the 
       xlab = "NDVI",
       xlim = c(0.2, 0.8),
       ylim = c(1, 60))
  dev.off()
}


##### Analyses #####
lm <- lm(FL016 ~ treatment, data = ___)
summary(lm)

lm <- lm(FL020 ~ treatment, data = ndvi)
summary(lm)

ndvi$ndvi_diff <- (ndvi$FL020_ndvi - ndvi$FL016_ndvi)
lm <- lm(ndvi_diff ~ treatment, data = ndvi)
summary(lm)






