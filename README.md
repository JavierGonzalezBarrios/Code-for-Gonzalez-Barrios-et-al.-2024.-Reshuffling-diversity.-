
#-------------------------------------------------------------------------------------------------------
## Gonzalez-Barrios et al., 2024. Nature Communications
## Emergent patterns of reef fish diversity correlate with coral assemblage shifts along the Great Barrier Reef
#-------------------------------------------------------------------------------------------------------

library(dplyr) 
library(tidyr)
library(tidyverse) 
library(stringr) 
library(ggplot2)
library(ggpubr)
library(marginaleffects)
library(ggeffects)
library(mgcv)
library(mgcv)
library(lme4)
library(betapart)

######################################################################
#-------------------- CHANGE IN SPECIES RICHNESS --------------------#

# set colors
mycolors <-c ("#fecc5c","#fd8d3c","#f03b20","#e31a1c","#67000d") 

# use multiple threads for fitting model
ctrl <- gam.control(nthreads = 4)

# model 1
m1 <- gam(s_mean_yr ~ REEF_LAT + s(REEF_LAT, by=ranges_5, k=4) + 
           s(SHELF, bs = 're') + s(REEF_NAME, bs = 're'), method="REML", 
         correlation = corAR1(form = ~ year), data= years_range, control = ctrl)

# plots global pattern in the GBR
sr1 <- plot_predictions(m1, condition = list("REEF_LAT", "ranges_5"), vcov = T) +
  coord_flip() + 
  scale_color_manual(values = mycolors_2, name = "Time interval") +
  scale_fill_manual(values = mycolors_2, name = "Time interval") +
  xlab("") + ylab("Species richness") + theme_classic() +
  theme(legend.position = c(0.8, 0.2)) 

# plot comparisons among periods
s1 <- plot_comparisons(m1,
                      variables = list(ranges_5 = c("1995-2000", "2001-2005")), condition = "REEF_LAT") + # Define comparison btw "2006-2015", "2016-2022"
  coord_flip() + ylim(-7, 11) +
  ggtitle ("2001-2005 - 1995-2000") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + theme(plot.title = element_text(size=10)) +
  xlab("Latitude") + ylab("")

s2 <- plot_comparisons(m1,
                      variables = list(ranges_5 = c("1995-2000", "2006-2010")), condition = "REEF_LAT") + # Define comparison btw "2006-2015", "2016-2022"
  coord_flip() + ylim(-7, 11) + 
  ggtitle ("2006-2010 - 1995-2000") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + theme(plot.title = element_text(size=10)) +
  xlab("") + ylab("")

s3 <- plot_comparisons(m1,
                      variables = list(ranges_5 = c("1995-2000", "2011-2015")), condition = "REEF_LAT") + # Define comparison btw "2006-2015", "2016-2022"
  coord_flip() + ylim(-7, 11) + 
  ggtitle ("2011-2015 - 1995-2000") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + theme(plot.title = element_text(size=10)) +
  xlab("Latitude") + ylab("Comparison")

s4 <- plot_comparisons(m1,
                      variables = list(ranges_5 = c("1995-2000", "2016-2022")), condition = "REEF_LAT") + # Define comparison btw "2006-2015", "2016-2022"
  coord_flip() + ylim(-7, 11) + 
  ggtitle ("2016-2022 - 1995-2000") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + theme(plot.title = element_text(size=10)) +
  xlab("") + ylab("Comparison")

# ------------------------Map Fig 1a
#-----
install.packages("remotes")
remotes::install_github("https://github.com/open-AIMS/gisaimsr")

library(raster)
library(tidyverse)
library(ggspatial)
library(sf)
library(dataaimsr)
library(gisaimsr)
library(ggrepel)

#---import latitudinal sector layer
lat_sec <- read_sf("Latitudinal sectors shape file/ltms.shp")

lat_sec  <- lat_sec %>% st_set_crs(7844)

#---data GBR shape
gbr_feat # data
gbr_feat_2 <- gbr_feat # duplicate data to plot different layers in the same map
sf::st_crs(gbr_feat_2)
#-----
gbr_feat %>%
  dplyr::filter(FEAT_NAME != "Mainland") %>%
  ggplot(data = .) +
  geom_sf()

#------------------
no_reef <- c("Cay", "Rock", "Island", "Sand")

#--------------- coordinates of latitudinal sectors
sec_data <- data.frame(sec_name=c("Cooktown/ \nLizard Island", "Cairns", "Innisfail", "Townsville", 
                                  "Whitsunday", "Pompey", "Swain", "Capricorn\nBunker"),
                       lat = c(-15,   -16.2,   -17.2,   -18.1, -19.4, -20.15, -20.6, -22.6),
                       lon = c(146.6, 146.5, 147.2, 148.2, 150.7, 151.7,  152.45, 151.5))

#     GBR map
map <- gbr_feat %>%
  filter(!FEAT_NAME %in% no_reef) %>%
  ggplot(data = .) +
  geom_sf(mapping = aes(fill = FEAT_NAME), color = "grey30", fill= "grey50") +
  geom_sf(data = subset(gbr_feat, FEAT_NAME %in% c("Reef")), 
          mapping = aes(fill = FEAT_NAME), color = "grey", fill= "grey") +
  
  geom_text(data = sec_data, aes(x = lon, y = lat, label = sec_name), # add sectors names
            size = 3.9, col = "black") + 
  
  geom_sf(data = lat_sec, alpha= 0) + # latitudinal sector borders
  
  theme_classic() +
  labs(x= "Longitude", y= "Latitude", color="") +
  
  coord_sf(xlim = c(145, 152.9), ylim = c(-14, -24.5), expand = F) +
  
  annotation_scale(location = "bl", width_hint = 0.4, style = "ticks") +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_minimal) + 
  
  geom_point(data = data_new_1, aes(x = long, y = lat),
             size = 3, shape = 21, fill= "#56B4E9", color = "black") +
  
  scale_colour_manual("", values = c("gbr_feat"="grey50", "data_new_1"="#56B4E9"))

map 

#-------------------------------------
#Australia country map
aus <- ggplot(data = world) +
  geom_sf(color = "grey30", fill = "grey50") +
  coord_sf(xlim = c(114, 154), ylim = c(-11, -43), expand = T, datum = NA) +
  geom_rect(xmin = 145, xmax = 153, ymin = -13.5, ymax = -24.5, 
            fill = NA, colour = "black", linewidth = 0.5) +
  theme_void()
# aus

fi1_a <- map + annotation_custom(
  grob = ggplotGrob(aus),
  xmin = 148.5,
  xmax = 153.1,
  ymin = -14,
  ymax = -17.5) 

fi1_a
#-----------------------------------------


######################################################################
#-------------------- TRENDS IN SPECIES TURNOVER --------------------#

ex_1 <- data2 %>%
  dplyr::select(A_SECTOR, REEF_NAME, SHELF, curp, cREPORT_YEAR, FISH_CODE, Average_abun) %>% 
  pivot_wider(names_from = FISH_CODE, values_from = Average_abun) %>% 
  arrange(curp, cREPORT_YEAR) %>%        
  as_tibble()

ex_1[is.na(ex_1)] <- 0

length(unique(ex_1$curp)) #  306 `curp;`102 when 'reef name'

col1 <- c(5:201) # specify columns to change list type to numeric

ex_1[ , col1] <- apply(ex_1[ , col1], 2,            # Specify own function within apply
                       function(x) as.numeric(as.character(x)))
# sapply(ex_1, class)  # check data class


#-----------------------------------------------------------------------------------
#            reef fish beta diversity
#----------------------------------------------

# prior transform abundance data

ex_1_sqrt <- ex_1 %>%
  mutate_at(vars(5:201), sqrt)


# -------------- loop
# prior need to define empty matrix to beta outputs
base2 <- data.frame(REEF_NAME = ex_1_sqrt$REEF_NAME,
                          curp = ex_1_sqrt$curp, 
                          Year = ex_1_sqrt$cREPORT_YEAR,
                          SHELF = ex_1_sqrt$SHELF,
                          subs.beta.bray.bal=NA, 
                          subs.beta.bray.gra=NA, 
                          subs.beta.bray=NA,
                          succ.beta.bray.bal=NA, 
                          succ.beta.bray.gra=NA, 
                          succ.beta.bray=NA)

# Create loop vector for every Site
for(i in unique(ex_1_sqrt$curp)) {
  
  # Get row index for first value of each group e.g. "Year 1"
  x <- min(which(ex_1_sqrt$curp == i)) 
  # Counter, only needed for calculating "successive years"
  y <- 0
  
  # Create loop vector one less than site observations 
  # e.g. exclude "Year 1"
  for(j in 1:(length(which(ex_1_sqrt$curp == i))-1)) {
    
    beta.subs <- beta.pair.abund(ex_1_sqrt[c(x,x+j),5:201], index.family="bray")

    beta.matrix[x+j,4:6] <- do.call(cbind, beta.subs)
    
    beta.succ <- beta.pair.abund(ex_1_sqrt[c(x+y, x+j),5:201], index.family="bray")
    # Add result to beta.matrix
    beta.matrix[x+j,7:9] <- do.call(cbind, beta.succ)
    # Increase "counter" by 1 for each inner loop iteration
    y <- y + 1
    
  }
  
}

#set colors
mycolors_8<-c("#a50026", "#d73027", "#f46d43", "#fdae61", "#abd9e9", "#74add1", "#4575b4", "#313695")

# GAM for Year-to-Year turnover GBR level
m2 <- gam(succ.beta.bray ~ year + s(year, k=8) + s(SHELF, bs = 're') + s(REEF_NAME, bs = 're'),
          data= base2, method = 'REML', correlation = corAR1(form = ~ year))    

t1<-plot_predictions(m2, condition = list("year"), vcov = T) +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") + 
  ggtitle("Great Barrier Reef") + xlab("") + ylab("Reef fish turnover (β)") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_classic()

# GLMM for Year-to-Year turnover GBR level
m6 <- lme4::lmer(succ.beta.bray ~ ranges_5 + (1|REEF_NAME) + (1|SHELF), data= base2)

t2<- plot_predictions(m6, condition = list("ranges_5"), vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("") + 
  theme_classic() + theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#GAM for Year-to-Year turnover at sector level
m3 <- gam(succ.beta.bray ~ year + A_SECTOR + s(year, by = A_SECTOR, k=8, m=1) + 
            s(SHELF, bs = 're') + s(REEF_NAME, bs = 're'), data= base2, method = 'REML', control = ctrl)

t3 <- plot_predictions(m3, condition = list("year", "A_SECTOR"), vcov = T) +
  scale_color_manual(values = mycolors_8) +
  scale_fill_manual(values = mycolors_8, name="") +
  labs(color = "", shape = "Latitudinal\nsector", 
       linetype = "Latitudinal\nsector") +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("Latitudinal sector") + xlab("") + ylab("Reef fish turnover (β)")  +
  theme(plot.margin=unit(c(-0.5,.1,.1,.1), "cm")) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(0.15, 0.38)) +
  theme_classic() + theme(legend.position = "none")

t4<- plot_predictions(m3, condition = list("A_SECTOR")) +
  theme(legend.position = "none") +
  ggtitle("") + xlab("") + ylab("") +
  coord_flip() +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL")) +
  theme(plot.margin=unit(c(-0.5,.1,.1,.1), "cm")) +
  theme_classic()

#-------------------------------------------------------------------------------
# GAM for Reference-Year turnover GBR level
m4 <- gam(subs.beta.bray ~ year + s(year, k=8) + s(SHELF, bs = 're') + s(REEF_NAME, bs = 're'),
          data= base2, method = 'REML', correlation = corAR1(form = ~ year))   

t5<-plot_predictions(m4, condition = list("year"), vcov = T) +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") + 
  ggtitle("Great Barrier Reef") + xlab("") + ylab("Reef fish turnover (β)") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_classic()

# GLMM for Reference-Year turnover GBR level
m7 <- lme4::lmer(subs.beta.bray ~ ranges_5 + (1|REEF_NAME) + (1|SHELF), data= bas2)

t6<- plot_predictions(m7, condition = list("ranges_5"), vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("") + 
  theme_classic() + theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#GAM for Year-to-Year turnover at sector level
m5 <- gam(subs.beta.bray ~ year + A_SECTOR + s(year, by = A_SECTOR, k=8, m=1) + 
            s(SHELF, bs = 're') + s(REEF_NAME, bs = 're'), data= base, method = 'REML', control = ctrl)

t7 <- plot_predictions(m5, condition = list("year", "A_SECTOR"), vcov = T) +
  scale_color_manual(values = mycolors_8) +
  scale_fill_manual(values = mycolors_8, name="") +
  labs(color = "", shape = "Latitudinal\nsector", 
       linetype = "Latitudinal\nsector") +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("Latitudinal sector") + xlab("") + ylab("Reef fish turnover (β)")  +
  theme(plot.margin=unit(c(-0.5,.1,.1,.1), "cm")) +
  scale_y_continuous(labels = label_number(accuracy = 0.01), limits = c(0.15, 0.38)) +
  theme_classic() + theme(legend.position = "none")

t8 <- plot_predictions(m5, condition = list("A_SECTOR")) +
  theme(legend.position = "none") +
  ggtitle("") + xlab("") + ylab("") +
  coord_flip() +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL")) +
  theme(plot.margin=unit(c(-0.5,.1,.1,.1), "cm")) +
  theme_classic()

######################################################################
#------------------------ CORAL COVER TRENDS ------------------------#

# Trends in coral cover in the GBR
m8 <- gam(ave_cover ~ cREPORT_YEAR  + s(cREPORT_YEAR, k=8) + 
             s(REEF_NAME, bs = 're') + s(SHELF, bs = 're'),
           correlation = corAR1(form = ~ cREPORT_YEAR), 
           data= data3, method = 'REML', control = ctrl)

t9 <- plot_predictions(m8, condition= c("cREPORT_YEAR"), vcov = T, draw=F) + 
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("") + xlab("") + ylab("Coral cover (%)") + 
  scale_x_continuous(breaks = scales::breaks_width(5)) +
  theme_classic() 

# GLMM for coral cover and time ranges
m10 <-  lme4::lmer(ave_cover ~ ranges + (1|REEF_NAME) + (1|SHELF), data= data3)

t10 <- plot_predictions(m10, condition= c("ranges"), vcov = T) +
  ggtitle("") + xlab("") + ylab("") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



######################################################################
#---------------- SHIFTS IN CORAL COMPOSITION TRENDS ----------------#

# -------------- loop
# prior need to define empty matrix to beta outputs
data3 <- data.frame(REEF_NAME = ex_1_sqrt$REEF_NAME,
                    curp = ex_1_sqrt$curp, 
                    Year = ex_1_sqrt$cREPORT_YEAR,
                    SHELF = ex_1_sqrt$SHELF,
                    c.subs.beta.bray.bal=NA, 
                    c.subs.beta.bray.gra=NA, 
                    c.subs.beta.bray=NA,
                    c.succ.beta.bray.bal=NA, 
                    c.succ.beta.bray.gra=NA, 
                    c.succ.beta.bray=NA)

# Create loop vector for every Site
for(i in unique(ex_1_sqrt$curp)) {
  
  # Get row index for first value of each group e.g. "Year 1"
  x <- min(which(ex_1_sqrt$curp == i)) 
  # Counter, only needed for calculating "successive years"
  y <- 0
  
  # Create loop vector one less than site observations 
  # e.g. exclude "Year 1"
  for(j in 1:(length(which(ex_1_sqrt$curp == i))-1)) {
    
    beta.subs <- beta.pair.abund(ex_1_sqrt[c(x,x+j),5:45], index.family="bray")
    
    beta.matrix[x+j,4:6] <- do.call(cbind, beta.subs)
    
    beta.succ <- beta.pair.abund(ex_1_sqrt[c(x+y, x+j),5:45], index.family="bray")
    # Add result to beta.matrix
    beta.matrix[x+j,7:9] <- do.call(cbind, beta.succ)
    # Increase "counter" by 1 for each inner loop iteration
    y <- y + 1
    
  }
  
}


# Trends in shift in coral composition in the GBR
m9 <-  gam(c.succ.beta.bray ~ year + s(year, k=8) + 
             s(REEF_NAME, bs = 're') + s(SHELF, bs = 're'), control = ctrl,
           data= data3, method = 'REML', correlation = corAR1(form = ~ year)) 

t11 <- plot_predictions(m9, condition= c("cREPORT_YEAR"), vcov = T, draw=F) + 
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("") + xlab("") + ylab("Coral cover (%)") + 
  scale_x_continuous(breaks = scales::breaks_width(5)) +
  theme_classic() 

# GLMM for shift in coral composition and time intervals
m11 <-  lme4::lmer(c.succ.beta.bray ~ ranges_5 + 
                     (1|REEF_NAME) + (1|SHELF), data= data3)

t12 <- plot_predictions(m11, condition = list("ranges_5"), vcov = T) +
  ggtitle("") + xlab("") + ylab("") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 



##########################################################################
# FISH TURNOVER CORRELATE WITH CHANGE IN CORAL COVER AND SHIFT IN CORAL COMPOSITION #

m12 <- lme4::lmer(succ.beta.bray ~ abs_abs_rate_scale + c.succ.beta.bray_scale + 
                    (1|REEF_NAME) + (1|SHELF), data= data3)

t13 <- ggpredict(m12, c("abs_abs_rate_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Change in coral cover (%)") + ylab("Reef fish turnover (β)") + ggtitle("")

t14 <- ggpredict(m12, c("c.succ.beta.bray_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Shift in coral composition") + ylab("Reef fish turnover (β)") + ggtitle("")

t15 <- plot_model(m12) + ylim(0, 0.017) + theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("") + ylab("Effect size") + ggtitle("")  +
  scale_x_discrete(limit = c("c.succ.beta.bray_scale", "abs_abs_rate_scale"),
                   labels = c( "Shift in coral\ncomposition", "Change in\ncoral cover\n(%)"))

# Shift in coral composition by latitudinal sector
m14 <- lmer(succ.beta.bray ~ c.succ.beta.bray * A_SECTOR + (1|REEF_NAME) + (1|SHELF), data= data3)

t16 <- plot_slopes(m14, variables = "c.succ.beta.bray", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  xlab("") + ylab("Effect size") + ggtitle("Reef fish turnover (β)") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t17 <- plot_predictions(m14, newdata= data3, by=c("c.succ.beta.bray", "A_SECTOR"), vcov = T, points= 0.2) +
  ggtitle("") + xlab("Shift in coral composition") + ylab("Reef fish turnover (β)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() + facet_wrap(~ A_SECTOR, ncol=2) +
  theme(legend.position = "none")

# Reef fish species turnover by change in coral cover at latitudinal sector
m15 <- lmer(succ.beta.bray ~ abs_abs_rate * A_SECTOR + (1|REEF_NAME) + (1|SHELF), data= data3)

t18 <- plot_slopes(m15, variables = "abs_abs_rate", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  ggtitle("") + xlab("") + ylab("Effect sizes") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t19 <- plot_predictions(m15, newdata= data3, by=c("abs_abs_rate", "A_SECTOR"), vcov = T, points= 0.3) +
  ggtitle("") + xlab("Change in coral cover (%)") + ylab("Reef fish turnover (β)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() +
  facet_wrap(~A_SECTOR, ncol=2) +
  theme(legend.position = "none")

##################################################################################
# FISH SPECIES RICHNESS BY CHANGE IN CORAL COVER AND SHIFTS IN CORAL COMPOSITION #

m13 <- lmer(sr_rate ~ abs_abs_rate + c.succ.beta.bray + (1|REEF_NAME) + (1|SHELF), data= data3)

t20 <- ggpredict(m13, c("abs_abs_rate_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Change in coral cover (%)") + ylab("Δ Species richness (α)") + ggtitle("")

t21 <- ggpredict(m13, c("c.succ.beta.bray_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Shift in coral composition") + ylab("Δ Species richness (α)") + ggtitle("")

t22 <- plot_model(m13) + theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("") + ylab("Standarized effect size") + ggtitle("")  +
  scale_x_discrete(limit = c("c.succ.beta.bray_scale", "abs_abs_rate_scale"),
                   labels = c( "Shift in coral\ncomposition", "Change in\ncoral cover\n(%)"))

# Change in species richness by shift in coral composition at latitudinal sector
m16 <- lmer(sr_rate ~ c.succ.beta.bray * A_SECTOR + (1|REEF_NAME) + (1|SHELF), data= data3)

t23 <- plot_slopes(m16, variables = "c.succ.beta.bray", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  xlab("") + ylab("Effect sizes") + ggtitle("") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t24 <- plot_predictions(m16, newdata= fish_coral_2, by=c("c.succ.beta.bray", "A_SECTOR"), vcov = T, points= 0.3) +
  ggtitle("") + xlab("Shift in coral composition") + ylab("Δ Species richness (α)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() + facet_wrap(~ A_SECTOR, ncol=2) +
  theme(legend.position = "none")
