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
  geom_sf(data = plots_feature, size = 3, color = "black", fill = "cyan1") + 
  ggtitle("Plot Boundaries") + 
  coord_sf()






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

barplot(cbind(FL016, FL020) ~  treatment.x, data = all_mean)


## add error bars on barplot
data_sum <- ddply(all_mean, c("treatment.x"), summarise,
               N    = length(FL016),
               mean = mean(FL016),
               sd   = sd(FL016),
               se   = sd / sqrt(N)
) 

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
  


aov <- aov(FL016 ~  FL020 + treatment.x, data = all_mean)
summary(aov)

aov <- aov(FL020 ~ treatment.x, data = all_mean)
summary(aov)


##### subset for all treatments ####
all_mean$ndvi_diff <- all_mean$FL016 - all_mean$FL020
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
aov_ct <- aov(FL016 ~  FL020, data = CT)
summary(aov_ct)

aov_gr <- aov(FL016 ~  FL020, data = GR)
summary(aov_gr)

aov_sh <- aov(FL016 ~  FL020, data = SH)
summary(aov_sh)

aov_gs <- aov(FL016 ~  FL020, data = GS)
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



boxplot(FL016 ~ plot_num, data = all_new) 


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






