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
FL016 <- raster("georef_rasters/FL016_georef.tif")  
FL016

FL020 <- raster("georef_rasters/FL020_goeref.tif")
FL020

FL020_crop <- crop(FL020, FL016)
FL020b<- resample(FL020_crop, FL016)
diff <- FL020b - FL016

plot(diff, 
     xlab = "longitude", 
     ylab = "latitude", 
     main = "NDVI difference")

RGB <- stack("RU_CYN_TR1_FL016_RGB.tif")
RGB

plot(FL016)
plot(FL020)

GPS <- na.omit(read.csv("CYN_plot_centers.csv"))
GPS2 <- subset(GPS, select= -c(plot, elevation))

## re order columns so longitude is first
GPS_order <- GPS2[,c("longitude", "latitude")]

## turn into spatial points
GPS_order2 <- SpatialPoints(GPS_order,
                            proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
CRS.new <- CRS("+proj=aea +lat_1=50 +lat_2=70 +lat_0=56 +lon_0=100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
GPS_final<- spTransform(GPS_order2, CRS.new)


plot(FL016, 
     main = "Pre-Clipping NDVI", 
     xlab = "Longitude", 
     ylab = "Latitude")
points(GPS_final$longitude, GPS_final$latitude, pch=21,
       bg = "black", col = "yellow")


plot(FL020, 
     main = "Post-Clipping NDVI", 
     xlab = "Longitude", 
     ylab = "Latitude")
points(GPS_final$longitude, GPS_final$latitude, pch=21,
       bg = "black", col = "yellow")



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

plots_feature <- st_read("plots/plots.shp")


st_geometry_type(plots_feature)

plots_feature

ggplot() + 
  geom_sf(data = plots_feature, size = 1, color = "black", fill = "cyan1") + 
  ggtitle("Plot Boundaries") + 
  coord_sf()


plot(RGB)

plotRGB(RGB, r = 1, g = 2, b = 3)
plot(plots_feature, add = TRUE, zoom = 4)
zm()



###### extract MEAN ndvi by feature class ######
pre_mean<- extract(FL016, plots_feature,
                      fun = mean,
                      df = TRUE)
names(pre_mean)[names(pre_mean) == "FL016_georef"] <- "FL016"


post_mean<- extract(FL020, plots_feature,
                   fun = mean,
                   df = TRUE)
names(post_mean)[names(post_mean) == "FL020_goeref"] <- "FL020"


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
plot(FL016 ~  treatment.x, 
     data = all_mean, 
     xlab = "treatment", 
     ylab = "NDVI", 
     main = "Pre-clipping NDVI (per treatment)")
plot(FL020 ~ treatment.x, 
     data = all_mean,
     xlab = "treatment", 
     ylab = "NDVI", 
     main = "Post-clipping NDVI (per treatment)")



##### add error bars on barplot #####
install.packages("plyr")
library(plyr)


install.packages("doBy")
library(doBy)

sum <- summaryBy(FL016 + FL020 ~ treatment.x, 
                   data=all_mean, 
                   FUN=c(length,mean,sd))

ggplot(sum, aes(x=treatment.x, y=FL016.mean, colour=FL016.mean)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=FL016.mean-FL016.sd, ymax=FL016.mean+FL016.sd), width=.1) +
  geom_line() +
  geom_point()+
  xlab("treatment") +
  ylab("NDVI")+
  ggtitle("Pre-clipping NDVI")

ggplot(sum, aes(x=treatment.x, y=FL020.mean, colour=FL020.mean)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=FL020.mean-FL020.sd, ymax=FL020.mean+FL020.sd), width=.1) +
  geom_line() +
  geom_point() + 
  xlab("treatment") +
  ylab("NDVI")+
  ggtitle("Post-clipping NDVI")



##### creating bar plots side by side #####
library(reshape2)
all_mean2 <- subset(all_mean, select= -c(ID.y, ID.x, treatment.y, plot_num))
df <- melt(all_mean2, id.vars='treatment.x')

cdata <- ddply(df, c("treatment.x", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
               )

####### BARGRAPH #####
ggplot(cdata, aes(x=treatment.x, y=mean, fill=variable)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  xlab("treatment") +
  ylab("NDVI")+
  ggtitle("NDVI pre- and post- clipping") +
  scale_fill_hue(name="Variable", # Legend label, use darker colors
                 breaks=c("FL016", "FL020"),
                 labels=c("Pre-Clipping", "Post-Clipping"))
  
##### BOXPLOT #####
ggplot(df, aes(x=treatment.x, y=value, fill=variable)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("NDVI")+
  scale_x_discrete(labels= c("Control", "Grass", "Grass + Shrub", "Shrub")) + 
  ggtitle("NDVI pre- and post- clipping") +
  scale_fill_manual(values=c("#CC6666", "#9999CC"),
                 name="Variable", # Legend label, use darker colors
                 breaks=c("FL016", "FL020"),
                 labels=c("Pre-Clipping", "Post-Clipping")) 


aov <- aov(value ~  variable + treatment.x, data = df)
summary(aov)
aov




aov <- aov(FL016 ~ treatment.x, data = all_mean) ##not different
summary(aov)
TukeyHSD(aov)

aov <- aov(FL020 ~ treatment.x, data = all_mean) ## some different
summary(aov)
TukeyHSD(aov)


# linear model




##### subset for all treatments ####
all_mean$ndvi_diff <- all_mean$FL020 - all_mean$FL016
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


##### plot of NDVI difference across treatments #####
barplot(ndvi_diff ~ plot_num, 
        CT)
barplot(ndvi_diff ~ plot_num, 
        GR) 
barplot(ndvi_diff ~ plot_num, 
        SH) 
barplot(ndvi_diff ~ plot_num, 
        GS) 


##### ANOVAs #####
aov_ct <- aov(FL016 ~  FL020, data = CT) ## no sign diff (p=0.406)
summary(aov_ct)

aov_gr <- aov(FL016 ~  FL020, data = GR) ## sif diff!!! (p=0.0245)
summary(aov_gr)

aov_sh <- aov(FL016 ~  FL020, data = SH) ## no sig diff (p=0.129)
summary(aov_sh)

aov_gs <- aov(FL016 ~  FL020, data = GS) ## no sig diff (p=0.181)
summary(aov_gs)





##### correlation plots #####
install.packages("ggcorrplot")
install.packages("GGally")
install.packages("tidyr")

library(tidyr)
library(ggcorrplot)
library(GGally)




veg = all_pt_int %>% 
  dplyr::select(plot_num:percent_composition) %>% 
  pivot_wider(names_from = functional_groups, 
              values_from = percent_composition,
              values_fill = 0)

veg

names(veg)


veg_var = veg %>%
  dplyr::select(plot_num, FL016, CON:FORB)
veg_var

corr <- round(cor(veg_var), 1)
corr

ggcorrplot(corr)

veg_var2 = veg %>%
  dplyr::select(plot_num, FL020, CON:FORB)
veg_var2

corr2 <- round(cor(veg_var2), 1)
corr2

ggcorrplot(corr2)


##### extract ALL ndvi by feature class ####
detach("package:tidyr", unload = TRUE) # extract wont ork with tidyr
library(sf)

plots_feature <- st_read("plots/plots.shp")


st_geometry_type(plots_feature)

plots_feature


pre <- extract(FL016, plots_feature,
               fun = NULL,
                   small = FALSE,
                   df = TRUE)
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




########### biomass removal vs ndvi difference #######

bio_removal <- na.omit(read.csv("CYN_bio_removal.csv"))
names(bio_removal)[names(bio_removal) == "Plot.ID"] <- "plot_num"
bio_removal2 <- aggregate(bio_removal$BIOMASS..g., by = list(plot=bio_removal$plot_num), FUN = sum)
names(bio_removal2)[names(bio_removal2) == "x"] <- "bio_removed"
names(bio_removal2)[names(bio_removal2) == "plot"] <- "plot_num"

ndvi_br <- na.omit(merge(all_mean, bio_removal2, by = c("plot_num"), all.x = TRUE, all.y = TRUE)) 

barplot(bio_removed ~ plot_num, data = ndvi_br,
        main = "Biomass Removed by Plot",
        xlab = "Plot", 
        ylab = "Biomass Removed (g)", 
        ylim = c(0,600), 
        col = "lightgreen")
lm <- lm(bio_removed ~ ndvi_diff, data = ndvi_br)
summary(lm)


aov <- aov(bio_removed ~ ndvi_diff + treatment.x, data = ndvi_br)
summary(aov)

ttest <- t.test(ndvi_br$ndvi_diff, ndvi_br$treatment.x)
summary(aov)

##### create box plot for ndvi diff #####
install.packages("ggpubr")
library(ggpubr)


cdata2 <- ddply(ndvi_br, c("treatment.x"), summarise,
               N    = length(ndvi_diff),
               mean = mean(ndvi_diff),
               sd   = sd(ndvi_diff),
               se   = sd / sqrt(N)
)


ggplot(ndvi_br, aes(x=treatment.x, y=ndvi_diff, fill=treatment.x)) + 
  geom_boxplot() +
  scale_x_discrete(labels= c("Grass", "Grass + Shrub", "Shrub")) + 
  xlab("Treatment") +
  ylab("NDVI Difference")+
  ggtitle("") +
  scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"),
                    name="Treatmeents", # Legend label, use darker colors
                    breaks=c(),
                    labels=c())
   

aov <- aov(ndvi_diff ~ treatment.x, data = ndvi_br)
summary(aov)
TukeyHSD(aov)



##### scatter plot (need this?) #####
plot(bio_removed ~ ndvi_diff, data = ndvi_br, 
     xlab = "NDVI Difference", 
     ylab = "Biomass Removed (g)", 
     main = "Biomass Removed vs NDVI",
     pch = 19, 
     col = "red")
abline(lm)

##### biomass removal by treatment ######
GR <- subset(ndvi_br, treatment.x == "GR")
SH <- subset(ndvi_br, treatment.x == "SH")
GS <- subset(ndvi_br, treatment.x == "GS")

lm<- lm(bio_removed ~ ndvi_diff, data = GR)
summary(lm)


lm<- lm(bio_removed ~ ndvi_diff, data = SH)
summary(lm)


lm<- lm(bio_removed ~ ndvi_diff, data = GS)
summary(lm)

## plot of biomass removed vs ndvi by treatmend (colored) ##
plot(bio_removed ~ ndvi_diff, 
     data = GR,
     xlim = c(-0.2, 0.2),
     ylim = c(0, 500),
     xlab = "NDVI Difference", 
     ylab = "Biomass Removed (g)", 
     main = "NDVI Difference vs Biomass Removed", 
     pch = 19, 
     col = "green")
points(bio_removed ~ ndvi_diff, 
       data = SH, 
       pch = 19, 
       col = "red")
points(bio_removed ~ ndvi_diff, 
       data = GS, 
       pch = 19, 
       col = "purple")
legend("topleft", c("grass removed", "shrub removed", "grass & shrub removed"),
       fill = c("green", "red", "purple"))

#####  subset by functional type #####
con <- subset(all_pt_int, functional_groups == "CON")
evsh <- subset(all_pt_int, functional_groups == "EVSH")
desh <- subset(all_pt_int, functional_groups == "DESH")
gram <- subset(all_pt_int, functional_groups == "GRAM")
forb <- subset(all_pt_int, functional_groups == "FORB")
cwd <- subset(all_pt_int, functional_groups == "CWD")
moss <- subset(all_pt_int, functional_groups == "MOSS")
lichen <- subset(all_pt_int, functional_groups == "LICH")
brg <- subset(all_pt_int, functional_groups == "BRG")
litr <- subset(all_pt_int, functional_groups == "LITR")
equ <- subset(all_pt_int, functional_groups == "EQU")


lm <- lm(FL016 ~ percent_composition, data = con)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = con,
     main = "NDVI vs Percent Cover (Conifers)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)

lm <- lm(FL016 ~ percent_composition, data = evsh)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = evsh,
     main = "NDVI vs Percent Cover (Evergreen Shrub)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)


lm <- lm(FL016 ~ percent_composition, data = desh)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = desh,
     main = "NDVI vs Percent Cover (Deciduous Shrub)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)



lm <- lm(FL016 ~ percent_composition, data = gram)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = gram,
     main = "NDVI vs Percent Cover (Graminoid)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)



lm <- lm(FL016~ percent_composition, data = forb)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = forb,
     main = "NDVI vs Percent Cover (Forb)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)



lm <- lm(FL016 ~ percent_composition, data = cwd)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = cwd,
     main = "NDVI vs Percent Cover (CWD)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)



lm <- lm(FL016 ~ percent_composition, data = lichen)
summary(lm)

plot(FL016 ~ percent_composition, 
     data = lichen,
     main = "NDVI vs Percent Cover (Lichen)",
     xlab = "Percent Cover", 
     ylab = "NDVI (pre-clipping)",
     pch = 19, 
     col = "darkgreen")
abline(lm)



## for loop to create histogram of NDVI values for each plot 

install.packages("reshape") 
library(reshape) 



for (i in 1:length(FL016_new)) { # for every column in the "new" data frame
  x <- FL016_new[ ,i] # identifying columns (?)
  # Plot histogram of x
  jpeg(file = paste("dist", names((FL016_new)[i]), ".jpeg", sep = ""))
  hist(x,
       main = paste("NDVI distribution for P", names((FL016_new)[i])), #paste name of column to the 
       xlab = "NDVI",
       xlim = c(0.2, 0.8),
       ylim = c(1, 110))
  dev.off()
}




##### NORMALIZING NDVI MAPS #####
diff_extract <- extract(diff, )

diff

cellStats(diff, stat = 'mean')
# [1] 0.03505925

norm <- FL020 - 0.03505925
plot(norm)
plot(FL016)

pre_mean2<- extract(FL016, plots_feature,
                   fun = mean,
                   df = TRUE)
names(pre_mean2)[names(pre_mean2) == "FL016_georef"] <- "FL016"


post_mean2<- extract(norm, plots_feature,
                    fun = mean,
                    df = TRUE)
names(post_mean2)[names(post_mean2) == "FL020_goeref"] <- "FL020"

plots_feature <- as.data.frame(plots_feature)
plots_feature$ID <- c(1:23)
plots = plots_feature %>% select(plot_num, ID)


pre_mean2 <- merge (pre_mean2, plots, by = "ID")
post_mean2 <- merge (post_mean2, plots, by = "ID")


treatments <- read.csv("treatments.csv")
treatments 
treatments$plot = gsub("P", "", treatments$plot)
names(treatments)[names(treatments) == "plot"] <- "plot_num"


pre_mean2 <- merge(pre_mean2, treatments, by = "plot_num")
post_mean2 <- merge(post_mean2, treatments, by = "plot_num")

all_mean_b <- merge(pre_mean2, post_mean2, by = "plot_num")
all_pt_int2 <-merge(all_mean_b, pt_int, by = "plot_num")

all_pt_int2$ndvi_diff <- all_pt_int2$FL016 - all_pt_int2$FL020

###### box plot #####
library(reshape2)
all_mean_b <- subset(all_mean_b, select= -c(ID.y, ID.x, treatment.y, plot_num))
df2 <- melt(all_mean_b, id.vars='treatment.x')

cdata2 <- ddply(df2, c("treatment.x", "variable"), summarise,
               N    = length(value),
               mean = mean(value),
               sd   = sd(value),
               se   = sd / sqrt(N)
)


##### BOXPLOT #####
ggplot(df2, aes(x=treatment.x, y=value, fill=variable)) + 
  geom_boxplot() +
  xlab("Treatment") +
  ylab("NDVI")+
  scale_x_discrete(labels= c("Control", "Grass", "Grass + Shrub", "Shrub")) + 
  ggtitle("NDVI pre- and post- clipping") +
  scale_fill_manual(values=c("#CC6666", "#9999CC"),
                    name="Variable", # Legend label, use darker colors
                    breaks=c("FL016", "FL020"),
                    labels=c("Pre-Clipping", "Post-Clipping")) 


aov <- aov(value ~  variable + treatment.x, data = df)
summary(aov)
aov

install.packages("nlme")
library(nlme)

lme = lme(value ~ variable+treatment.x, data = df2)



