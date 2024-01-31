


## Gonzalez-Barrios et al., 2024. 
## Multidecadal coral assemblage shifts have reshuffled reef fish diversity along a latitudinal gradient
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

######################################################################
#-------------------- CHANGE IN SPECIES RICHNESS --------------------#

# set colors
mycolors <-c ("#fecc5c","#fd8d3c","#f03b20","#e31a1c","#67000d") 

# use multiple threads for fitting model
ctrl <- gam.control(nthreads = 4)

# model 1
m1 <- gam(s_mean_yr ~ REEF_LAT + s(REEF_LAT, by=ranges_5, k=3), method="REML", 
         correlation = corAR1(form = ~ year), data= data1, control = ctrl)

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


######################################################################
#-------------------- TRENDS IN SPECIES TURNOVER --------------------#

#set colors
mycolors_8<-c("#a50026", "#d73027", "#f46d43", "#fdae61", "#abd9e9", "#74add1", "#4575b4", "#313695")

# GAM for Year-to-Year turnover GBR level
m2 <- gam(succ.beta.bray ~ year + s(year, k=8) + s(REEF_NAME, bs = 're'),
                data= base2, method = 'REML', correlation = corAR1(form = ~ year))   

t1<-plot_predictions(m2, condition = list("year"), vcov = T) +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") + 
  ggtitle("Great Barrier Reef") + xlab("") + ylab("Reef fish turnover (β)") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_classic()

# GLMM for Year-to-Year turnover GBR level
m6 <- lme4::lmer(succ.beta.bray ~ ranges_5 + (1|REEF_NAME), data= base2)

t2<- plot_predictions(m6, condition = list("ranges_5"), vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("") + 
  theme_classic() + theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#GAM for Year-to-Year turnover at sector level
m3 <- gam(succ.beta.bray ~ year + A_SECTOR + s(year, by = A_SECTOR, k=8, m=1) + 
            s(REEF_NAME, bs = 're'), data= base2, method = 'REML', control = ctrl)

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
m4 <- gam(sub.beta.bray ~ year + s(year, k=8) + s(REEF_NAME, bs = 're'),
          data= base2, method = 'REML', correlation = corAR1(form = ~ year))   

t5<-plot_predictions(m4, condition = list("year"), vcov = T) +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") + 
  ggtitle("Great Barrier Reef") + xlab("") + ylab("Reef fish turnover (β)") +
  scale_y_continuous(labels = label_number(accuracy = 0.01)) +
  theme_classic()

# GLMM for Year-to-Year turnover GBR level
m7 <- lme4::lmer(sub.beta.bray ~ ranges_5 + (1|REEF_NAME), data= base2)

t6<- plot_predictions(m7, condition = list("ranges_5"), vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("") + 
  theme_classic() + theme(legend.position = "none",
                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

#GAM for Year-to-Year turnover at sector level
m5 <- gam(sub.beta.bray ~ year + A_SECTOR + s(year, by = A_SECTOR, k=8, m=1) + 
            s(REEF_NAME, bs = 're'), data= base2, method = 'REML', control = ctrl)

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

# GLMM for coral cover and time ranges
m8 <-  lme4::lmer(ave_cover ~ ranges + (1|REEF_NAME), data= data3)

t9 <- plot_predictions(m8, condition= c("ranges"),  vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("")

# Trends in coral cover in the GBR
t10 <- data3 %>%
  ggplot(aes(cREPORT_YEAR, ave_cover)) + 
  stat_summary(geom = "line", fun.y = mean, size=1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("") + xlab("") + ylab("Coral cover (%)") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020))

######################################################################
#---------------- SHIFTS IN CORAL COMPOSITION TRENDS ----------------#

# GLMM for shift in coral composition and time ranges
m9 <-  lme4::lmer(c.succ.beta.bray ~ ranges_5 + (1|REEF_NAME), data= data3)

t11 <- plot_predictions(m9, condition = list("ranges_5"), vcov = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("") + xlab("") + ylab("")

# Trends in shift in coral composition in the GBR
t12 <-  data3 %>%
  ggplot(aes(year, c.succ.beta.bray)) + 
  stat_summary(geom = "line", fun.y = mean, size=1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.3) +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  geom_vline(xintercept=c(2000, 2005, 2010, 2015), linetype="dashed") +
  ggtitle("") + xlab("") + ylab("Shift in coral composition") + 
  scale_x_continuous(breaks = c(1995, 2000, 2005, 2010, 2015, 2020))

##########################################################################
# FISH TURNOVER BY CHANGE IN CORAL COVER AND SHIFTS IN CORAL COMPOSITION #

m10 <- lme4::lmer(succ.beta.bray ~ abs_abs_rate_scale + c.succ.beta.bray_scale + (1|REEF_NAME),
                  data= data3)

t13 <- ggpredict(m10, c("abs_abs_rate_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Change in coral cover (%)") + ylab("Reef fish turnover (β)") + ggtitle("")

t14 <- ggpredict(m10, c("c.succ.beta.bray_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Shift in coral composition") + ylab("Reef fish turnover (β)") + ggtitle("")

t15 <- plot_model(m10) + ylim(0, 0.017) + theme_classic() + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("") + ylab("Standarized effect size") + ggtitle("")  +
  scale_x_discrete(limit = c("c.succ.beta.bray_scale", "abs_abs_rate_scale"),
                   labels = c( "Shift in coral\ncomposition", "Change in\ncoral cover\n(%)"))

# Shift in coral composition by latitudinal sector
m12 <- lmer(succ.beta.bray ~ c.succ.beta.bray * A_SECTOR + (1|REEF_NAME), data= data3)

t16 <- plot_slopes(m12, variables = "c.succ.beta.bray", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  xlab("") + ylab("Effect size") + ggtitle("Reef fish turnover (β)") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t17c<- plot_predictions(m12, newdata= fish_coral_2, by=c("c.succ.beta.bray", "A_SECTOR"), vcov = T, points= 0.2) +
  ggtitle("") + xlab("Shift in coral composition") + ylab("Reef fish turnover (β)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() + facet_wrap(~ A_SECTOR, ncol=2) +
  theme(legend.position = "none")

# Reef fish species turnover by change in coral cover at latitudinal sector
m13 <- lmer(succ.beta.bray ~ abs_abs_rate * A_SECTOR + (1|REEF_NAME), data= data3)

t18c<- plot_slopes(tl4, variables = "abs_abs_rate", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  ggtitle("") + xlab("") + ylab("Effect sizes") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t19 <- plot_predictions(tl4, newdata= fish_coral_2, by=c("abs_abs_rate", "A_SECTOR"), vcov = T, points= 0.3) +
  ggtitle("") + xlab("Change in coral cover (%)") + ylab("Reef fish turnover (β)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() +
  facet_wrap(~A_SECTOR, ncol=2) +
  theme(legend.position = "none")

##################################################################################
# FISH SPECIES RICHNESS BY CHANGE IN CORAL COVER AND SHIFTS IN CORAL COMPOSITION #

m11 <- lmer(sr_rate ~ abs_abs_rate + c.succ.beta.bray + (1|REEF_NAME), data= data3)

t20 <- ggpredict(m11, c("abs_abs_rate_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Change in coral cover (%)") + ylab("Δ Species richness (α)") + ggtitle("")

t21 <- ggpredict(m11, c("c.succ.beta.bray_scale")) %>%
  plot(add.data=T, line.size=1, dot.size=1.5) + theme_classic() + 
  xlab("Shift in coral composition") + ylab("Δ Species richness (α)") + ggtitle("")

t22 <- plot_model(m11) + theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("") + ylab("Standarized effect size") + ggtitle("")  +
  scale_x_discrete(limit = c("c.succ.beta.bray_scale", "abs_abs_rate_scale"),
                   labels = c( "Shift in coral\ncomposition", "Change in\ncoral cover\n(%)"))

# Change in species richness by shift in coral composition at latitudinal sector
m14 <- lmer(sr_rate ~ c.succ.beta.bray * A_SECTOR + (1|REEF_NAME), data= fish_coral_2)

t23 <- plot_slopes(m14, variables = "c.succ.beta.bray", condition ="A_SECTOR")  + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme_classic() + coord_flip() +
  xlab("") + ylab("Effect sizes") + ggtitle("") +
  scale_x_discrete(limits=c("CB", "SW", "PO", "WH", "TO","IN", "CA", "CL"))

t24 <- plot_predictions(m14, newdata= fish_coral_2, by=c("c.succ.beta.bray", "A_SECTOR"), vcov = T, points= 0.3) +
  ggtitle("") + xlab("Shift in coral composition") + ylab("Δ Species richness (α)") +
  scale_color_manual(values = mycolors_8, name="") +
  scale_fill_manual(values = mycolors_8, name="") +
  theme_classic() + facet_wrap(~ A_SECTOR, ncol=2) +
  theme(legend.position = "none")
