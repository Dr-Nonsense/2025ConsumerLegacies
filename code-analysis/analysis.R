# created 120.03.2025 by Seraina Cappelli
# rerun CWM analyses 2009-2024 biomass based calculation, linear models only

# Bonaparte, Sophie Hunger - Daft Punk Spielen in meinem Haus


## SETUP ####

library(tidyverse)
# dplyr     1.1.4     # readr     2.1.5
# forcats   1.0.0     # stringr   1.5.1
# ggplot2   3.5.1     # tibble    3.2.1
# lubridate 1.9.4     # tidyr     1.3.1
# purrr     1.0.2  
library(broom)      # 1.14.0
library(lme4)       # 1.1-35.1
library(lmerTest)   # 3.1-3
library(MuMIn)      # 1.47.5
library(kableExtra) # 1.4.0
library(cowplot)    # 1.1.2
library(viridis)    # 0.6.4
library(ggeffects)  # 1.3.4
library(emmeans)    # 1.10.5
library(ggpubr)     # 0.6.0
library(ggrepel)    # 0.9.6
library(vegan)      # 2.6-4
library(ggtext)     # 0.1.2
library(glue)       # 1.8.0

# project path:
EnemyRemovalDIR <- getwd() %>% paste(., "/", sep = "")


# prep some things to be used in figures
theme_set(theme_cowplot()+
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))+
            theme(axis.text = element_text(size=8),
                  axis.title = element_text(size=10),
                  legend.text = element_text(size=8),
                  legend.title = element_text(size=10),
                  title = element_text(size=10),
                  text = element_text(size = 8)))

myvalues <- c("grey", "#3399FF",  "#44AA99", "#FF9999", "#882255", "black")  
mybreaks <- c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides")
mylabels <- c("Control", "Fenced", "Insecticide", "Soil drench fungicide", "Foliar fungicide", "All pesticides")

traitlabels <- list(Nmass        = "CWM Tissue N", 
                    DiasporeMass = "CWM Seed mass",
                    PlantHeight  = "CWM Plant height",
                    LeafArea     = "CWM Leaf area", 
                    LMA          = "CMW LMA",
                    LDMC         = "CWM LDMC",
                    invsimpson   = "Inverse Simpson",
                    PC1          = "PC1",
                    PC2          = "PC2")

traitlabels_unit <- list(Nmass   = "CWM Tissue N [%]", 
                         DiasporeMass = "CWM Seed mass [mg]",
                         PlantHeight  = "CWM Plant height [m]",
                         LeafArea     = "CWM Leaf area ["~mm^2~"]", 
                         LMA          = "CMW LMA [g "~m^-2~"]",
                         LDMC         = "CWM LDMC",
                         invsimpson   = "Inverse Simpson",
                         PC1          = "PC1",
                         PC2          = "PC2")
traitlabels_unit_2 <- list(Nmass        = "Treatment effect on CWM Tissue N [%]", 
                           DiasporeMass = "Treatment effect on CWM Seed mass [mg]",
                           PlantHeight  = "Treatment effect on CWM Plant height [m]",
                           LeafArea     = "Treatment effect on CWM Leaf area ["~mm^2~"]", 
                           LMA          = "Treatment effect on CMW LMA [g "~m^-2~"]",
                           LDMC         = "Treatment effect on CWM LDMC",
                           invsimpson   = "Treatment effect on Inverse Simpson",
                           PC1          = "PC1",
                           PC2          = "PC2")

ybreaks <- seq(1, 16, 5)
ylabels <- seq(2009, 2024, 5)

sigbreaks <- c("p > 0.05", "p < 0.05")
sigvalues <- c("dashed", "solid")
sigvalues2 <- c(0.5, 1)

# set up contrasts to use in the models
# Biodiversity Experiment
BBc_Control.ER       <- matrix(c(1, -0.25, -0.25, -0.25, -0.25)) # control vs. enemies removed
BBc_Single.All       <- matrix(c(0,   1/3,   1/3,  1/3,    -1)) # single vs. all enemies removed
BBc_Herbivores.Fungi <- matrix(c(0,     1,  -0.5, -0.5,     0)) # herbivores vs. fungi
BBc_Soil.Foliar      <- matrix(c(0,     0,     1,   -1,     0)) # soil vs. foliar (within fungi)

BBc <- matrix(c(BBc_Control.ER, BBc_Single.All, BBc_Herbivores.Fungi, BBc_Soil.Foliar), ncol = 4)
colnames(BBc) <-c(": Control vs. ER", ": Single vs. All", ": Herbivores vs. Fungi", ": Soil vs. Foliar")

# Oldfield Experiment
OFc_Control.ER       <- matrix(c(1, -0.2, -0.2, -0.2, -0.2, -0.2)) # control vs. enemies removed
OFc_Single.All       <- matrix(c(0, 0.25, 0.25, 0.25, 0.25, -1))   # single vs. all enemies removed
OFc_Herbivores.Fungi <- matrix(c(0,  0.5,  0.5, -0.5, -0.5,  0))   # herbivores vs. fungi
OFc_Insects.Mammals  <- matrix(c(0,    1,   -1,    0,    0,  0))   # insects vs. mammal (within herbivores)
OFc_Soil.Foliar      <- matrix(c(0,    0,    0,    1,   -1,  0))   # soil vs. foliar (within fungi)

OFc_Large.Small      <- matrix(c(0,    1, -1/3, -1/3, -1/3,  0))   # large vs. small (instead of herbivores vs. fungi)
OFc_Insect.Fungi     <- matrix(c(0,    0,    1, -0.5, -0.5,  0))   # Insects vs. fungi (instead of insects vs. mammal)


OFc1 <- matrix(  c(OFc_Control.ER,         OFc_Single.All,     OFc_Herbivores.Fungi,     OFc_Insects.Mammals,    OFc_Soil.Foliar), ncol = 5)
colnames(OFc1) <-c(": Control vs. ER", ": Single vs. All", ": Herbivores vs. Fungi", ": Insects vs. Mammals", ": Soil vs. Foliar")

OFc2 <- matrix(  c(OFc_Control.ER,         OFc_Single.All,    OFc_Large.Small,       OFc_Insect.Fungi,    OFc_Soil.Foliar), ncol = 5)
colnames(OFc2) <-c(": Control vs. ER", ": Single vs. All", ": Large vs. Small", ": Insects vs. Fungi", ": Soil vs. Foliar")

rm(BBc_Control.ER, BBc_Single.All, BBc_Herbivores.Fungi, BBc_Soil.Foliar, 
   OFc_Control.ER, OFc_Single.All, OFc_Herbivores.Fungi, OFc_Insects.Mammals, 
   OFc_Soil.Foliar, OFc_Large.Small, OFc_Insect.Fungi)

## GET DATA ####
CWM_Bio <- read.csv(paste(EnemyRemovalDIR, "/data-derived/e244_CWMTraits.csv", sep = "")) %>%
  filter(!Year %in% c(2007)) %>%
  mutate(logNumSp = log(NumSp),
         Year1 = Year-2008,
         Treatment = factor(Treatment,
                            levels = c("Control", 
                                       "Insecticide", 
                                       "SoilDrenchFungicide", 
                                       "FoliarFungicide",
                                       "AllPesticides"))) %>%
  
  merge(read.csv(paste(EnemyRemovalDIR, "/data-derived/E244_div.csv", sep = "")) %>%
          select(Plot, Year, invsimpson_biomass) %>%
          rename("CWM_invsimpson_Biom_sownSp" = "invsimpson_biomass") %>%
          filter(!Year %in% c(2007)),
        by = c("Plot", "Year"),
        all = T) %>%
  
  merge(read.csv(paste(EnemyRemovalDIR, "/data-derived/e244_CWMFGAbundance.csv", sep = "")) %>%
          select(Year, Plot, Functional.group, Mass.g.m.2.) %>%
          filter(!Year %in% c(2007)) %>%
          filter(Functional.group %in% c("L", "F", "C3", "C4")) %>%
          mutate(Mass.g.m.2. = replace_na(Mass.g.m.2., 0)) %>%
          pivot_wider(id_cols = c(Year, Plot),
                      values_from = Mass.g.m.2.,
                      values_fill = 0,
                      names_from = Functional.group) %>%
          rename("L_biom" = "L",
                 "F_biom" = "F",
                 "C3_biom" = "C3",
                 "C4_biom" = "C4"),
        by = c("Year", "Plot"), all = T) %>%
  
  filter(!NumSp %in% 1) %>%
  filter(!Year %in% 2008)


CWM_OF <- read.csv(paste(EnemyRemovalDIR, "/data-derived/e245_CWMTraits.csv", sep = ""))  %>%
  filter(!Year %in% c(2007)) %>%
  filter(!Subplot %in% "West") %>%
  mutate(Year1 = Year-2008,
         Treatment = factor(Treatment,
                            levels = c("Control", 
                                       "Fenced",
                                       "Insecticide", 
                                       "SoilDrenchFungicide", 
                                       "FoliarFungicide",
                                       "AllPesticides")))%>%
  
  merge(read.csv(paste(EnemyRemovalDIR, "/data-derived/E245_div.csv", sep = "")) %>%
          filter(!Subplot %in% "West") %>%
          select(Plot, Year, invsimpson_biomass) %>%
          rename("CWM_invsimpson_Biom" = "invsimpson_biomass") %>%
          filter(!Year %in% c(2007)),
        by = c("Plot", "Year"),
        all = T)%>%
  
  merge(read.csv(paste(EnemyRemovalDIR, "/data-derived/e245_CWMFGAbundance.csv", sep = "")) %>%
          filter(!Year %in% c(2007)) %>%
          filter(!Subplot %in% "West") %>%
          select(Year, Plot, Functional.group, Mass.g.m.2.) %>%
          filter(Functional.group %in% c("L", "F", "C3", "C4")) %>%
          mutate(Mass.g.m.2. = replace_na(Mass.g.m.2., 0)) %>%
          pivot_wider(id_cols = c(Year, Plot),
                      values_from = Mass.g.m.2.,
                      values_fill = 0,
                      names_from = Functional.group) %>%
          rename("L_biom" = "L",
                 "F_biom" = "F",
                 "C3_biom" = "C3",
                 "C4_biom" = "C4"),
        by = c("Year", "Plot"), all = T) %>%
  filter(!Year %in% 2008) %>%
  mutate(Block = ceiling(Plot/6))


sp_change <- read.csv(paste(EnemyRemovalDIR, "/data-derived/species_biomass_change.csv", sep = "")) %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", 
                                                  "Fenced", 
                                                  "SoilDrenchFungicide", 
                                                  "Insecticide", 
                                                  "FoliarFungicide",
                                                  "AllPesticides")))
sp_TrtEff <- read.csv(paste(EnemyRemovalDIR, "/data-derived/species_biomass_TrtEff.csv", sep = ""))%>%
  mutate(Treatment = factor(Treatment, levels = c("Control", 
                                                  "Fenced", 
                                                  "SoilDrenchFungicide", 
                                                  "Insecticide", 
                                                  "FoliarFungicide",
                                                  "AllPesticides"))) %>%
  rename("LeafArea" = "Leaf_Area",
         "PlantHeight" = "Plant_height",
         "DiasporeMass" = "Diaspore_mass") %>%
  filter(!contrast %in% NA)


#___________________________________________________________________________####

## PCA & PERMANOVA ####
### PCA ####
CWM_pca <- 
  prcomp( 
    rbind.data.frame(
      CWM_Bio %>% 
        filter(!Year %in% c(2019, 2020)) %>%
        select(ends_with(c("_Biom_sownSp"))) %>%
        select(!CWM_invsimpson_Biom_sownSp)  %>%
        na.omit() %>% 
        rename_with(~str_remove(., '_Biom_sownSp')),
      CWM_OF %>%
        filter(!Year %in% c(2019, 2020)) %>%
        filter(!Subplot %in% "West") %>%
        select(ends_with("_Biom"))  %>%
        select(starts_with("CWM_")) %>%
        select(!CWM_invsimpson_Biom) %>%
        na.omit() %>% 
        rename_with(~str_remove(., '_Biom'))
    ),
    scale = T
  ) 


CWM_both <- rbind.data.frame(
  CWM_Bio %>%
    filter(!Year %in% c(2019, 2020)) %>%
    mutate(experiment = "Biodiversity") %>%
    select(experiment, Plot, Treatment, Year, NumSp, ends_with(c("_Biom_sownSp"))) %>%
    select(!CWM_invsimpson_Biom_sownSp)  %>%
    rename_with(~str_remove(., '_Biom_sownSp')) %>%
    na.omit(),
  CWM_OF %>%
    filter(!Year %in% c(2019, 2020)) %>%
    filter(!Subplot %in% "West") %>%
    mutate(experiment = "Oldfield",
           NumSp = "nope") %>%
    select(experiment, Plot, Treatment, Year, NumSp, starts_with("CWM_")) %>%
    select(experiment, Plot, Treatment, Year, NumSp, ends_with("_Biom")) %>%
    select(!CWM_invsimpson_Biom) %>%
    rename_with(~str_remove(., '_Biom')) %>%
    na.omit()
) %>%
  # introduces NAs because Oldfield doesnt have species richness manipulated
  mutate(NumSp = as.numeric(paste(NumSp)))



plot_data1 <- CWM_pca %>%
  augment(CWM_both) %>% # add original dataset back in
  rename(PC1 = .fittedPC1,
         PC2 = .fittedPC2,
         PC3 = .fittedPC3,
         PC4 = .fittedPC4)

plot_data2 <- CWM_pca %>% 
  tidy(matrix = "rotation") %>%  
  filter(PC %in% c(1,2,3,4)) %>% 
  pivot_wider(names_from = "PC", values_from = value, names_prefix = "PC") %>% 
  mutate(PC1 = PC1 * 7, PC2 = PC2 * 7, PC3 = PC3 * 7, PC4 = PC4 * 7)  %>%
  mutate(column = gsub(column, pattern = "CWM_", replacement = "")) %>%
  mutate(column = case_when(column %in% "LeafArea" ~ "Leaf area",
                            column %in% "PlantHeight" ~ "Plant\nheight",
                            column %in% "DiasporeMass" ~ "Seed mass",
                            column %in% "LMA" ~ "LMA",
                            column %in% "LDMC" ~ "LDMC",
                            column %in% "Nmass" ~ "Tissue N"))


##### -> FIGURE 1a ####
compplot <- ggplot(plot_data1, aes(x=PC1, y=PC2)) + 
  theme(legend.position = "bottom") +
  labs(x = paste("PC1 (", 100*round((CWM_pca %>% summary())$importance[2,1], 2), "%)", sep = ""),
       y = paste("PC2 (", 100*round((CWM_pca %>% summary())$importance[2,2], 2), "%)", sep = ""),
       color = "") +
  scale_color_manual(breaks = c("Oldfield", "Biodiversity"),
                     values = c("#DDBB55",  "#9988CC"),
                     labels = c("Oldfield experiment", "Biodiversity experiment"))+ 
  geom_point(size = 1.5, aes(color = experiment), alpha = 0.5) +
  geom_segment(data = plot_data2, xend = 0, yend = 0, arrow = arrow(type = "closed", end = "first", length = unit(0.03, "npc")), size = 1) +
  geom_text(data = plot_data2 ,  aes(label = column), hjust = c(0, 1, 0, 0, 1, 0), nudge_x = c(+0.1, -0.1, +0.1, +0.1, -0.1, +0.1), fontface = "bold", # size = 2.108759
            size = 4)+ 
  stat_ellipse(type = "norm", aes(color =experiment), show.legend = FALSE, size = 1.5)


#### trajectory through time ####
centers2 = data.frame(PC1 = numeric(), 
                      PC2 = numeric(), 
                      PC1_prev = numeric(), 
                      PC2_prev = numeric(), 
                      Year = integer(), 
                      experiment = factor(),
                      Treatment = factor())
for(i in unique(plot_data1$experiment)){
  subs1 <- plot_data1 %>% filter(experiment %in% i)
  
  for(k in unique(subs1$Treatment)){
    subs2 <- subs1 %>% filter(Treatment %in% k)
    
    
    for(j in 1:length(sort(unique(subs2$Year)))){
      subs3 <- subs2 %>% filter(Year %in% (sort(unique(subs2$Year)))[j])
      subs4 <- subs2 %>% filter(Year %in% (sort(unique(subs2$Year)))[j-1])
      xx <- cov.wt(subs3[,c("PC1", "PC2")])$center %>% 
        t() %>% 
        data.frame()
      yy <-  cov.wt(subs4[,c("PC1", "PC2")])$center %>% 
        t() %>% 
        data.frame() %>%
        rename("PC1_prev" = "PC1",
               "PC2_prev" = "PC2")
      xx  <- cbind.data.frame(xx, yy, Year = (sort(unique(subs2$Year)))[j], experiment = i, Treatment = k)
      
      centers2 = rbind(centers2, xx)  }
  }
}
centers2 <- centers2 %>%
  mutate(PC1_prev = case_when(Year %in% 2009 ~ NA, .default = PC1_prev),
         PC2_prev = case_when(Year %in% 2009 ~ NA, .default = PC2_prev))

# use linear models for PC1 and PC2 to estimate start and end point of trajectory through time
centers_trajectory <- data.frame()
for (i in mybreaks){
  df_bio <- centers2 %>% filter(Treatment %in% i & experiment %in% "Biodiversity")
  
  if (nrow(df_bio) > 0) {
    lm_PC1_bio <- lm(PC1 ~ Year, df_bio)
    lm_PC2_bio <- lm(PC2 ~ Year, df_bio)
    
    df_new_bio <- expand_grid(Year = c(2009, 2024),
                              experiment = "Biodiversity",
                              Treatment = i)
    df_new_bio$PC1 <- predict(lm_PC1_bio,  newdata = df_new_bio)
    df_new_bio$PC2 <- predict(lm_PC2_bio,  newdata = df_new_bio) 
    
    centers_trajectory <- rbind(centers_trajectory, df_new_bio)
  }
  
  df_OF <- centers2 %>% filter(Treatment %in% i & experiment %in% "Oldfield")
  
  if (nrow(df_OF) > 0) {
    lm_PC1_OF <- lm(PC1 ~ Year, df_OF)
    lm_PC2_OF <- lm(PC2 ~ Year, df_OF)
    
    df_new_OF <-  expand_grid(Year = c(2009, 2024),
                              experiment = "Oldfield",
                              Treatment = i)
    df_new_OF$PC1 <- predict(lm_PC1_OF,  newdata = df_new_OF)
    df_new_OF$PC2 <- predict(lm_PC2_OF,  newdata = df_new_OF)
    
    centers_trajectory <- rbind(centers_trajectory, df_new_OF)
  }
}

#### -> FIGURE S1 ####
compplot_trajectory_details <-ggplot(plot_data1 %>%
                                       mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))), 
                                     aes(PC1, PC2)) + 
  facet_grid(experiment~ Treatment,
             labeller = labeller(Treatment = c("Control" = "Control",
                                               "Fenced" = "Fenced",
                                               "Insecticide" = "Insecticide",
                                               "SoilDrenchFungicide" = "Soil Drench Fungicide",
                                               "FoliarFungicide" = "Foliar Fungicide",
                                               "AllPesticides" = "All Pesticides"),
                                 experiment = c("Biodiversity" = "Biodiversity experiment",
                                                "Oldfield" = "Oldfield experiment"))) +
  theme(legend.position = "bottom") +
  labs(x = paste("PC1 (", 100*round((CWM_pca %>% summary())$importance[2,1], 2), "%)", sep = ""),
       y = paste("PC2 (", 100*round((CWM_pca %>% summary())$importance[2,2], 2), "%)", sep = ""),
       color = "") +
  scale_color_viridis(breaks = ylabels, direction = -1) +
  geom_segment(data = expand_grid(plot_data2, 
                                  Treatment = c("Control", "Fenced",
                                                "Insecticide", 
                                                "SoilDrenchFungicide", 
                                                "FoliarFungicide", "AllPesticides"),
                                  experiment = c("Biodiversity", "Oldfield")) %>%
                 mutate(PC1 = PC1/4,
                        PC2 = PC2/4)%>%
                 mutate(Treatment = factor(Treatment, 
                                           levels = c("Control", "Fenced", 
                                                      "Insecticide", 
                                                      "SoilDrenchFungicide", 
                                                      "FoliarFungicide", "AllPesticides"))) %>%
                 filter(!(Treatment %in% "Fenced" & experiment %in% "Biodiversity")),
               xend  = 0, 
               yend  = 0, 
               arrow = arrow(type   = "closed", 
                             end    = "first", 
                             length = unit(0.03, "npc")), 
               size  = 0.5, color = "grey") +
  geom_text(data = expand_grid(plot_data2, 
                               Treatment = c("Control", "Fenced", "Insecticide", 
                                             "SoilDrenchFungicide", "FoliarFungicide",
                                             "AllPesticides"),
                               experiment = c("Biodiversity", "Oldfield")) %>%
              mutate(PC1 = PC1/4,
                     PC2 = PC2/4,
                     hjuster = ifelse(column %in% c("Tissue N", "Seed mass"), 1, 0))%>%
              mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced",
                                                              "Insecticide", 
                                                              "SoilDrenchFungicide", 
                                                              "FoliarFungicide", 
                                                              "AllPesticides")))%>%
              filter(!(Treatment %in% "Fenced" & experiment %in% "Biodiversity")),
            aes(label = column, hjust = hjuster),
            nudge_x = c(rep(0.1, 11), rep(-0.1, 11), rep(0.1, 22), rep(-0.1, 11), rep(0.1,11)),
            # fontface = "bold",
            size = 2.108759, color = "grey") +
  # size = 4) +
  geom_point(data = centers2%>%
               mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))), 
             aes(y = PC2, x = PC1, color = Year)) + 
  # ggrepel::geom_text_repel(data = centers2%>%
  #                            mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))),
  # aes(y = PC2, x = PC1, color = Year, label = Year)) + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = centers2 %>% filter(!Year %in% 2008)%>%
                 mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))),
               aes(yend = PC2, xend = PC1, y = PC2_prev, x = PC1_prev, color = Year), 
               arrow=arrow(length=unit(0.2,"cm"),  type = "closed"), size = 0.6)


#### -> FIGURE 1b ####
compplot_trajectory_details3 <-  ggplot(plot_data1 %>% 
           mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))), 
       aes(PC1, PC2)) + 
  facet_grid(experiment~.,
             labeller = labeller(experiment = c("Biodiversity" = "Biodiversity experiment",
                                                "Oldfield" = "Oldfield experiment"))) +
  theme(legend.position = "bottom") +
  labs(x = paste("PC1 (", 100*round((CWM_pca %>% summary())$importance[2,1], 2), "%)", sep = ""),
       y = paste("PC2 (", 100*round((CWM_pca %>% summary())$importance[2,2], 2), "%)", sep = "")) +
  scale_color_manual(breaks=mybreaks, values=myvalues, labels=mylabels) +
  scale_alpha_continuous(breaks=ybreaks+2008) +
  geom_segment(data = expand_grid(plot_data2, 
                                  Treatment = c("Control", "Fenced",
                                                "Insecticide", 
                                                "SoilDrenchFungicide", 
                                                "FoliarFungicide", "AllPesticides"),
                                  experiment = c("Biodiversity", "Oldfield")) %>%
                 mutate(PC1 = PC1/4,
                        PC2 = PC2/4)%>%
                 mutate(Treatment = factor(Treatment, 
                                           levels = c("Control", "Fenced", 
                                                      "Insecticide", 
                                                      "SoilDrenchFungicide", 
                                                      "FoliarFungicide", "AllPesticides"))) %>%
                 filter(!(Treatment %in% "Fenced" & experiment %in% "Biodiversity")),
               xend  = 0, 
               yend  = 0, 
               arrow = arrow(type   = "closed", 
                             end    = "first", 
                             length = unit(0.03, "npc")), 
               size  = 0.5, color = "lightgrey") +
  geom_text(data = expand_grid(plot_data2, 
                               Treatment = c("Control", "Fenced", "Insecticide", 
                                             "SoilDrenchFungicide", "FoliarFungicide",
                                             "AllPesticides"),
                               experiment = c("Biodiversity", "Oldfield")) %>%
              mutate(PC1 = PC1/4,
                     PC2 = PC2/4,
                     hjuster = ifelse(column %in% c("Tissue N", "Seed mass"), 1, 0))%>%
              mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced",
                                                              "Insecticide", 
                                                              "SoilDrenchFungicide", 
                                                              "FoliarFungicide", 
                                                              "AllPesticides")))%>%
              filter(!(Treatment %in% "Fenced" & experiment %in% "Biodiversity")),
            aes(label = column, hjust = hjuster, size = 4),
            nudge_x = c(rep(0.1, 11), rep(-0.1, 11), rep(0.1, 22), rep(-0.1, 11), rep(0.1,11)),
            fontface = "bold",
            size = 4, color = "lightgrey") +
  geom_point(data = centers2%>%
               mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))), 
             aes(y = PC2, x = PC1, color = Treatment, alpha = Year), size = 2) + 
  # ggrepel::geom_text_repel(data = centers2%>%
  #                            mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))),
  # aes(y = PC2, x = PC1, color = Year, label = Year)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color ="lightgrey") +
  geom_vline(xintercept = 0, linetype = "dotted", color ="lightgrey") +
  geom_segment(data = centers_trajectory %>%
                 pivot_wider(id_cols = c(experiment, Treatment), values_from = c(PC1, PC2), names_from= Year)%>%
                 mutate(Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))),
               aes(yend = PC2_2024, xend = PC1_2024, y = PC2_2009, x = PC1_2009, color = Treatment),
               arrow=arrow(length=unit(0.25,"cm"),  type = "closed"), size = 1.5, show.legend = FALSE)


rm(plot_data2, i, subs1, subs2, centers2, xx, yy, j,
   subs3, subs4, k, centers_trajectory)


### PERMANOVA ####
# PERMANOVAS with scaled data & Euclidian 
#### Biodiversity experiment ####
CWM_bio_traits <- CWM_Bio %>% 
  filter(!Year %in% c(2019, 2020)) %>%
  select(ends_with("_Biom_sownSp"), BigBioPlot, Plot, Treatment, Year, Year1, NumSp) %>% 
  select(!CWM_invsimpson_Biom_sownSp)  %>%
  rename_with(~str_remove(., '_Biom_sownSp')) %>%
  na.omit() %>% 
  mutate(CWM_LeafArea     = scale(CWM_LeafArea),
         CWM_Nmass        = scale(CWM_Nmass), 
         CWM_LMA          = scale(CWM_LMA),
         CWM_PlantHeight  = scale(CWM_PlantHeight),
         CWM_DiasporeMass = scale(CWM_DiasporeMass),
         CWM_LDMC         = scale(CWM_LDMC)) %>%
  select(starts_with("CWM"))

CWM_bio_env <- CWM_Bio %>% 
  filter(!Year %in% c(2019, 2020)) %>%
  select(ends_with("_Biom_sownSp"), BigBioPlot, Plot, Treatment, Year, Year1, NumSp) %>% 
  select(!CWM_invsimpson_Biom_sownSp)  %>%
  rename_with(~str_remove(., '_Biom_sownSp')) %>%
  na.omit() %>% 
  mutate(CWM_LeafArea     = scale(CWM_LeafArea),
         CWM_Nmass        = scale(CWM_Nmass), 
         CWM_LMA          = scale(CWM_LMA),
         CWM_PlantHeight  = scale(CWM_PlantHeight),
         CWM_DiasporeMass = scale(CWM_DiasporeMass),
         CWM_LDMC         = scale(CWM_LDMC)) %>%
  select(BigBioPlot, Plot, Treatment, Year, Year1, NumSp) %>%
  mutate(logNumSp = log(NumSp)) %>%
  mutate(Year3 = as.factor(Year))

### -> Only plot as RANEF, no interactions
permanova_bio <- with(CWM_bio_env, 
                      adonis2(CWM_bio_traits ~ Treatment + Year1 + logNumSp, 
                              data = CWM_bio_env, 
                              permutations = 20, 
                              strata = Plot,
                              method = "euclidean")) %>%
  data.frame() %>%
  rename("p-value" = "Pr..F.") %>%
  rownames_to_column("Variable") %>%
  mutate(Variable = case_when(Variable %in% "Year1" ~ "Year",
                              Variable %in% "logNumSp" ~ "log species richness",
                              .default = Variable))


#### Oldfield experiment ####
CWM_of_traits <- CWM_OF %>%
  filter(!Year %in% c(2019, 2020)) %>%
  filter(!Subplot %in% "West") %>%
  select(ends_with("_Biom"), Plot, Treatment, Year) %>%
  select(starts_with("CWM_"), Plot, Treatment, Year) %>%
  select(!CWM_invsimpson_Biom) %>% 
  rename_with(~str_remove(., '_Biom')) %>%
  na.omit() %>%
  mutate(CWM_LeafArea     = scale(CWM_LeafArea),
         CWM_Nmass        = scale(CWM_Nmass), 
         CWM_LMA          = scale(CWM_LMA),
         CWM_PlantHeight  = scale(CWM_PlantHeight),
         CWM_DiasporeMass = scale(CWM_DiasporeMass),
         CWM_LDMC         = scale(CWM_LDMC)) %>%
  select(starts_with("CWM_"))
CWM_of_env <- CWM_OF %>%
  filter(!Year %in% c(2019, 2020)) %>%
  filter(!Subplot %in% "West") %>%
  select(ends_with("_Biom"), Plot, Treatment, Year, Year1) %>%
  select(starts_with("CWM_"), Plot, Treatment, Year, Year1) %>%
  select(!CWM_invsimpson_Biom) %>% 
  rename_with(~str_remove(., '_Biom')) %>%
  na.omit() %>%
  mutate(CWM_LeafArea     = scale(CWM_LeafArea),
         CWM_Nmass        = scale(CWM_Nmass), 
         CWM_LMA          = scale(CWM_LMA),
         CWM_PlantHeight  = scale(CWM_PlantHeight),
         CWM_DiasporeMass = scale(CWM_DiasporeMass),
         CWM_LDMC         = scale(CWM_LDMC)) %>%
  select(Plot, Treatment, Year, Year1) %>%
  mutate(Block = floor(Plot/6)*6,
         Year3 = factor(Year))

### --> Only plot as RANEF, no interactions
permanova_OF <- with(CWM_of_env, 
                     adonis2(CWM_of_traits ~ Treatment + Year1, 
                             data = CWM_of_env, 
                             permutations = 20, 
                             strata = Plot,
                             method = "euclidean")) %>%
  data.frame() %>%
  rename("p-value" = "Pr..F.") %>%
  rownames_to_column("Variable") %>%
  mutate(Variable = case_when(Variable %in% "Year1" ~ "Year",
                              .default = Variable))


rm(CWM_bio_traits, CWM_bio_env, CWM_of_traits, CWM_of_env)


#___________________________________________________________________________####

## MODELS Biodiversity Experiment ####
mods_Bio <- list()
FG_mods_Bio <- list()
anova_Bio <- data.frame(Var = c(),
                        trait = c(), 
                        X2    = c(),
                        DF    = c(),
                        p     = c(),
                        P_print   = c(),
                        X2_print  = c(),
                        stars     = c(),
                        Var_print = c())

anova_for_print_Bio <- list()
anova_for_print_FG_Bio <- list()
mod_summary_Bio <- list()
treatment_effects_Bio <- list()
plots_Bio <- list()
plots_FG_Bio <- list()


for (i in c("Nmass", "DiasporeMass", "PlantHeight", "LeafArea", "LMA", "LDMC", "invsimpson", "PC1", "PC2")){
  
  ### data ####
  df <- CWM_Bio %>%
    merge(.,
          plot_data1 %>% 
            filter(experiment %in% "Biodiversity") %>% 
            select(Plot, Treatment, Year, PC1, PC2) %>%
            rename("CWM_PC1_Biom_sownSp"= "PC1",
                   "CWM_PC2_Biom_sownSp" = "PC2"),
          by = c("Plot", "Treatment", "Year"), all.x = T) %>% 
    select(BigBioPlot, Plot, Treatment, Year, Year1, logNumSp, L_biom, F_biom, C3_biom, C4_biom, ends_with("sownSp")) %>%
    select(BigBioPlot, Plot, Treatment, Year, Year1, logNumSp, L_biom, F_biom, C3_biom, C4_biom, contains(i)) %>%
    rename("y_b" = paste("CWM", i, "Biom_sownSp", sep = "_"))  %>%
    mutate(Treatment=`contrasts<-`(factor(Treatment), , BBc)) %>%
    # inverse simpson is Inf in monoculture
    filter(!y_b %in% Inf)
  
  
  ### models ####
  lin_mod_b = lmer(y_b ~ Treatment * logNumSp * Year1 + (1|BigBioPlot/Plot) + (1|Year), data = df)
  mods_Bio[[i]] <- lin_mod_b
  
  lin_mod_b_L  <- update(lin_mod_b, .~. + L_biom)
  lin_mod_b_F  <- update(lin_mod_b, .~. + F_biom)
  lin_mod_b_C3 <- update(lin_mod_b, .~. + C3_biom)
  lin_mod_b_C4 <- update(lin_mod_b, .~. + C4_biom)
  
  FG_mods_Bio[[i]] <- list(base_model   = lin_mod_b,
                           legume_model = lin_mod_b_L,
                           forb_model   = lin_mod_b_F,
                           C3_model     = lin_mod_b_C3,
                           C4_model     = lin_mod_b_C4)
  
  ### anova table (df & for print) ####
  anova_Bio <- rbind.data.frame(anova_Bio,
                                car::Anova(lin_mod_b) %>% 
                                  data.frame() %>%
                                  rename("X2" = "Chisq",
                                         "DF" = "Df",
                                         "p" = "Pr..Chisq.") %>%
                                  rownames_to_column("Var") %>%
                                  mutate(trait = i)  %>%
                                  mutate(P_print2 = case_when(p < 0.001 ~ "p < 0.001", .default = paste("p =", round(p, digits = 3))),
                                         P_print = case_when(p < 0.001 ~ "<0.001", .default = paste(round(p, digits = 3))),
                                         X2_print = round(X2, digits = 4),
                                         stars = case_when(p < 0.001 ~ "***",
                                                           p > 0.001 & p <0.01 ~ "**",
                                                           p > 0.01  & p < 0.05 ~ "*",
                                                           p > 0.05 & p <0.1 ~ ".",
                                                           .default = ""),
                                         Var_print = case_when(Var %in% "Treatment" ~ "Treatment (T)",
                                                               Var %in% "logNumSp"  ~ "log Species Richness (SR)",
                                                               Var %in% "Year1"     ~ "Year (Y)",
                                                               Var %in% "I(Year1^2)" ~ "Year^2 (Y2)",
                                                               Var %in% "Treatment:logNumSp"   ~ "T x SR",
                                                               Var %in% "Treatment:Year1"      ~ "T x Y",
                                                               Var %in% "Treatment:I(Year1^2)" ~ "T x Y2",
                                                               Var %in% "logNumSp:Year1"       ~ "SR x Y",
                                                               Var %in% "logNumSp:I(Year1^2)"  ~ "SR x Y2",
                                                               Var %in% "Treatment:logNumSp:Year1"      ~ "T x SR x Y",
                                                               Var %in% "Treatment:logNumSp:I(Year1^2)" ~ "T x SR x Y2"),
                                         print = paste(Var_print, ": ", P_print2, stars, sep = "")) 
                                ) 
  
  
  anova_for_print_Bio[[i]] <- anova_Bio %>%
    filter(trait %in% i) %>%
    select(!trait) %>%
    mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
    mutate(P_print = gsub(P_print,
                          pattern = '<span style=" font-weight: bold; " >0</span>',
                          replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
    select(Var_print, X2_print, DF, P_print) %>%
    kbl(escape = F, 
        caption = paste("Biodiversity experiment:", traitlabels[[i]], sep = " ") , 
        col.names = c("Term", "X2", "DF",  "P-value"))  %>%
    kable_minimal(full_width = F, html_font = "Cambria") %>%
    kable_paper()
  
  
  #### -> Tables S3,5,7,9,11,13,15 ####
  AIC_header <-  rbind.data.frame("AIC", 
                                  AIC(lin_mod_b, lin_mod_b_L, lin_mod_b_F, lin_mod_b_C3, lin_mod_b_C4) %>%
                                              round(0) %>%
                                              data.frame %>%
                                              select(AIC)) %>%
    cbind.data.frame(c(1, rep(2, 5)))
  
  
  anova_for_print_FG_Bio[[i]] <- car::Anova(lin_mod_b) %>% 
    data.frame() %>%
    rownames_to_column("Var") %>%
    mutate(model = "base_model") %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_L) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "legume_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_F) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "forb_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_C3) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "C3_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_C4) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "C4_model")) %>%
    
    rename("X2" = "Chisq",
           "p" = "Pr..Chisq.") %>%
    mutate(Var = case_when(Var %in% c("L_biom", "F_biom", "C3_biom", "C4_biom") ~  "FG_biomass", .default = Var)) %>%
    mutate(P_print = case_when(p < 0.001 ~ "<0.001", .default = paste(round(p, digits = 3))),
           X2_print = round(X2, digits = 4),
           Var_print = case_when(Var %in% "Treatment" ~ "Treatment (T)",
                                 Var %in% "logNumSp"  ~ "log Species Richness (SR)",
                                 Var %in% "Year1"     ~ "Year (Y)",
                                 Var %in% "Treatment:logNumSp"   ~ "T x SR",
                                 Var %in% "Treatment:Year1"      ~ "T x Y",
                                 Var %in% "logNumSp:Year1"       ~ "SR x Y",
                                 Var %in% "Treatment:logNumSp:Year1"      ~ "T x SR x Y",
                                 Var %in% "FG_biomass" ~ "FG biomass")) %>%
    mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
    mutate(P_print = gsub(P_print,
                          pattern = '<span style=" font-weight: bold; " >0</span>',
                          replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
    select(Var_print, model, X2_print, P_print) %>%
 
    pivot_wider(id_cols = Var_print,
                values_from = c(X2_print, P_print),
                names_from = model) %>%
    select(Var_print, 
           X2_print_base_model,  P_print_base_model,
           X2_print_legume_model, P_print_legume_model,              
           X2_print_forb_model,  P_print_forb_model,
           X2_print_C3_model, P_print_C3_model,
           X2_print_C4_model, P_print_C4_model)  %>%
    
    kbl(escape = F, 
        caption = paste("Biodiversity experiment:", traitlabels[[i]], sep = " "),
        col.names = c("Term", rep(c("X2", "P-value"), 5)))  %>%
    add_header_above(AIC_header) %>%
    add_header_above(c(" ", "Base model" = 2, "Legume model" = 2, "Forb model" = 2, "C3 model" = 2, "C4 model" = 2))%>%
    kable_paper()
    
  
  
  
  ### model summary for print ####
  t1 <- lin_mod_b %>% 
    summary() %>%
    coef() %>% 
    data.frame() %>%
    round(digits = 4) %>%
    rename("p" = "Pr...t..") %>%
    mutate(p = cell_spec(p, bold = ifelse(p < 0.05, TRUE, FALSE))) %>%
    mutate(p = gsub(p, pattern = '<span style=" font-weight: bold; " >0</span>', replacement = '<span style=" font-weight: bold;    " ><0.001</span>')
    ) %>%
    rename("SE" = "Std..Error",
           "DF" = "df",
           "t-value" = "t.value",
           "p-value" = "p") %>%
    kbl(escape = F, 
        caption = paste("Biodiversity experiment: CWM", i, sep = " ")) %>%
    kable_minimal(full_width = F, html_font = "Cambria") %>%
    kable_paper()
  
  t2 <- VarCorr(lin_mod_b) %>%
    as.data.frame() %>%
    select(grp, vcov, sdcor) %>%
    mutate(vcov = round(vcov, digits=4),
           sdcor = round(sdcor, digits = 4)) %>%
    rename("Random effect" = "grp",
           "Variance" = "vcov",
           "DF" = "sdcor") %>%
    kbl(escape = F) %>%
    kable_minimal(full_width = F,  html_font = "Cambria") %>%
    kable_paper()
  
  mod_summary_Bio[[i]] <- list(coefs = t1,
                               ranefs = t2)
  
  ### emmeans pairwise ####
  # treatment_effects_Bio[[i]] <- (emmeans(lin_mod_b,  pairwise ~ Treatment))$contrasts %>%
  #   data.frame()
  
  ### plots ####
  #### -> FIGURE S2-4,9a ####
  fig_TY <- expand_grid(Year1 = c(1:16),
              Treatment = c("Control",
                            "Insecticide",
                            "SoilDrenchFungicide",
                            "FoliarFungicide",
                            "AllPesticides")) %>%
    mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
           Year = Year1+2008) %>%
    mutate(fixed_lb = predict(mods_Bio[[i]],  newdata = . , re.form = NA)) %>%

    merge(
      expand_grid(Year1 = c(1:16),
                  Treatment = c("Control",
                                "Insecticide",
                                "SoilDrenchFungicide",
                                "FoliarFungicide",
                                "AllPesticides")) %>%
        mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
               Year = Year1+2008) %>%
        filter(!Year %in% c(2019, 2020)) %>%
        mutate(full_lb = predict(mods_Bio[[i]],  newdata = . , re.form = ~(1|Year))),

      all = T, by = c("Year1", "Treatment", "logNumSp", "Year")
    ) %>%

    ggplot(aes(y = fixed_lb, x = Year, color = Treatment)) +
    labs(title = "Biodiversity experiment",
         y = traitlabels_unit[[i]] )+
    scale_color_manual(breaks = mybreaks, labels= mylabels, values= myvalues)+
    scale_x_continuous(breaks = ylabels) +
    geom_line() +
    geom_line(aes(y=full_lb), linetype="dotted") +
    geom_label(label = (anova_Bio %>%
                 filter(trait %in% i) %>%
                 filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
                 select(print) %>%
                 summarize(print = paste(print, collapse = "\n")))$print,
               color = "black",
               x = -Inf, hjust = 0,
               y = Inf, vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)
  
  #### -> FIGURE S2-4,9e ####
  fig_SRY <- ggpredict(mods_Bio[[i]], terms = c("Year1 [1:16, by = 5]", "logNumSp [1.3:3.6 by =0.2]")) %>%
    data.frame() %>%
    rename("Year1" = "x",
           "fixed_lb" = "predicted",
           "logNumSp" = "group") %>%
    mutate(logNumSp = as.numeric(paste(logNumSp)),
           NumSp = exp(logNumSp)) %>%
    select(Year1, NumSp, fixed_lb) %>%
    
    ggplot(aes(y = fixed_lb, x = NumSp, color = Year1, group = Year1)) +
    labs(title = "Biodiversity experiment",
         y = traitlabels_unit[[i]], 
         x = "Species richness") +
    scale_color_viridis(breaks = ybreaks,
                        labels = ylabels,
                        direction = -1)+
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    geom_line() +
    geom_label(label = (anova_Bio %>%
                          filter(trait %in% i) %>%
                          filter(Var %in% c("logNumSp", "Year1",  "logNumSp:Year1")) %>%
                          select(print) %>%
                          summarize(print = paste(print, collapse = "\n")))$print,
               color = "black",
               x = -Inf, hjust = 0,
               y = Inf, vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)
  
  #### -> FIGURE S2-4,9c ####
  fig_SRT <- ggpredict(mods_Bio[[i]], terms = c("Treatment", "logNumSp [1.3:3.6 by =0.2]")) %>%
    data.frame() %>%
    rename("Treatment" = "x",
           "fixed_lb" = "predicted",
           "logNumSp" = "group") %>%
    mutate(logNumSp = as.numeric(paste(logNumSp)),
           NumSp = exp(logNumSp)) %>%
    select(Treatment, NumSp, fixed_lb) %>%
    
    ggplot(aes(y = fixed_lb, x = NumSp, color = Treatment)) +
    labs(title = "Biodiversity experiment",
         y = traitlabels_unit[[i]], 
         x = "Species richness") +
    scale_color_manual(values=myvalues, breaks=mybreaks, labels=mylabels)+
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    geom_line() +
    geom_label(label = (anova_Bio %>%
                          filter(trait %in% i) %>%
                          filter(Var %in% c("logNumSp", "Treatment",  "Treatment:logNumSp")) %>%
                          select(print) %>%
                          summarize(print = paste(print, collapse = "\n")))$print,
               color = "black",
               x = -Inf, hjust = 0,
               y = Inf, vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)
  
  fig_TYSR <- ggpredict(mods_Bio[[i]], terms = c("logNumSp [1.3:3.6 by =0.2]", "Treatment", "Year1 [1:16, by = 5]")) %>% 
    data.frame() %>%
    mutate(facet = as.numeric(paste(facet)) + 2008) %>%
    ggplot(aes(y = predicted, x = exp(x), color = group, fill = group)) +
    theme(legend.position = "bottom")+
    facet_grid(.~facet) +
    scale_color_manual(breaks = mybreaks,
                       values = myvalues,
                       labels = mylabels) +
    scale_fill_manual(breaks = mybreaks,
                      values = myvalues,
                      labels = mylabels) +
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    labs(title = "Biodiversity experiment",
         y = traitlabels_unit[[i]],
         x = "Species richness",
         color = "Treatment",
         fill = "Treatment") +
    geom_line() +
    geom_label(data = cbind.data.frame(label = (anova_Bio %>% filter(trait %in% i & Var %in% "Treatment:logNumSp:Year1"))$print,
                                       x = -Inf,
                                       y = Inf,
                                       facet = 2009),
               aes(y = y, x = x, label = label),
               color = "black",
               # x = -Inf,
               hjust = 0,
               # y = Inf,
               vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)
  
  plots_Bio[[i]] <- list(Treatment_Year = fig_TY,
                         SR_Year = fig_SRY,
                         Treatment_SR = fig_SRT,
                         Treatment_SR_Year = fig_TYSR)
  
  
  
  
  ### plots of functional group model with lowest AIC ####

  
  # find functional group model with lowest AIC
  lowest_AIC_model <- cbind.data.frame(
    model = c("legume_model", "forb_model", "C3_model", "C4_model"),
    var = c("L_biom", "F_biom", "C3_biom", "C4_biom"),
    AIC =  (AIC(FG_mods_Bio[[i]]$legume_model,
                FG_mods_Bio[[i]]$forb_model,
                FG_mods_Bio[[i]]$C3_model,
                FG_mods_Bio[[i]]$C4_model))$AIC) %>%
    filter(AIC %in% min(AIC))
  
  
  
  ### plots
  #### -> FIGURE S2-4,9b ####
  fig_FG_TY <- expand_grid(Year1 = c(1:16),
                           Treatment = c("Control",
                                         "Insecticide",
                                         "SoilDrenchFungicide",
                                         "FoliarFungicide",
                                         "AllPesticides"),
                           !!lowest_AIC_model$var := mean(CWM_Bio[,lowest_AIC_model$var], na.rm =T)
  ) %>%
    mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
           Year = Year1+2008) %>%
    mutate(fixed_lb = predict(FG_mods_Bio[[i]][[lowest_AIC_model$model]],  newdata = . , re.form = NA)) %>%
    
    merge(
      expand_grid(Year1 = c(1:16),
                  Treatment = c("Control",
                                "Insecticide",
                                "SoilDrenchFungicide",
                                "FoliarFungicide",
                                "AllPesticides"),
                  !!lowest_AIC_model$var := mean(CWM_Bio[,lowest_AIC_model$var], na.rm =T)) %>%
        mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
               Year = Year1+2008) %>%
        filter(!Year %in% c(2019, 2020)) %>%
        mutate(full_lb = predict(FG_mods_Bio[[i]][[lowest_AIC_model$model]],  newdata = . , re.form = ~(1|Year))),
      
      all = T, by = c("Year1", "Treatment", "logNumSp", "Year")
    ) %>%
    
    ggplot(aes(y = fixed_lb, x = Year, color = Treatment)) +
    labs(title = paste("Biodiversity experiment,", gsub(lowest_AIC_model$model, pattern = "_", replacement = " ")),
         y = traitlabels_unit[[i]] )+
    scale_color_manual(breaks = mybreaks, labels= mylabels, values= myvalues)+
    scale_x_continuous(breaks = ylabels) +
    geom_line() +
    geom_line(aes(y=full_lb), linetype="dotted") #+
  # geom_label(label = (anova_Bio %>%
  #                       filter(trait %in% i) %>%
  #                       filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
  #                       select(print) %>%
  #                       summarize(print = paste(print, collapse = "\n")))$print,
  #            color = "black",
  #            x = -Inf, hjust = 0,
  #            y = Inf, vjust = 1,
  #            label.size = 0,
  #            fill = NA,
  #            size = 2.108759)
  
  #### -> FIGURE S2-4,9f ####
  fig_FG_SRY <- ggpredict(FG_mods_Bio[[i]][[lowest_AIC_model$model]], terms = c("Year1 [1:16, by = 5]", "logNumSp [1.3:3.6 by =0.2]")) %>%
    data.frame() %>%
    rename("Year1" = "x",
           "fixed_lb" = "predicted",
           "logNumSp" = "group") %>%
    mutate(logNumSp = as.numeric(paste(logNumSp)),
           NumSp = exp(logNumSp)) %>%
    select(Year1, NumSp, fixed_lb) %>%
    
    ggplot(aes(y = fixed_lb, x = NumSp, color = Year1, group = Year1)) +
    labs(title = paste("Biodiversity experiment,", gsub(lowest_AIC_model$model, pattern = "_", replacement = " ")),
         y = traitlabels_unit[[i]], 
         x = "Species richness") +
    scale_color_viridis(breaks = ybreaks,
                        labels = ylabels,
                        direction = -1)+
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    geom_line() #+
  # geom_label(label = (anova_Bio %>%
  #                       filter(trait %in% i) %>%
  #                       filter(Var %in% c("logNumSp", "Year1",  "logNumSp:Year1")) %>%
  #                       select(print) %>%
  #                       summarize(print = paste(print, collapse = "\n")))$print,
  #            color = "black",
  #            x = -Inf, hjust = 0,
  #            y = Inf, vjust = 1,
  #            label.size = 0,
  #            fill = NA,
  #            size = 2.108759)
  
  
  #### -> FIGURE S2-4,9d ####
  fig_FG_SRT <- ggpredict(FG_mods_Bio[[i]][[lowest_AIC_model$model]], terms = c("Treatment", "logNumSp [1.3:3.6 by =0.2]")) %>%
    data.frame() %>%
    rename("Treatment" = "x",
           "fixed_lb" = "predicted",
           "logNumSp" = "group") %>%
    mutate(logNumSp = as.numeric(paste(logNumSp)),
           NumSp = exp(logNumSp)) %>%
    select(Treatment, NumSp, fixed_lb) %>%
    
    ggplot(aes(y = fixed_lb, x = NumSp, color = Treatment)) +
    labs(title = paste("Biodiversity experiment,", gsub(lowest_AIC_model$model, pattern = "_", replacement = " ")),
         y = traitlabels_unit[[i]], 
         x = "Species richness") +
    scale_color_manual(values=myvalues, breaks=mybreaks, labels=mylabels)+
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    geom_line() #+
  # geom_label(label = (anova_Bio %>%
  #                       filter(trait %in% i) %>%
  #                       filter(Var %in% c("logNumSp", "Treatment",  "Treatment:logNumSp")) %>%
  #                       select(print) %>%
  #                       summarize(print = paste(print, collapse = "\n")))$print,
  #            color = "black",
  #            x = -Inf, hjust = 0,
  #            y = Inf, vjust = 1,
  #            label.size = 0,
  #            fill = NA,
  #            size = 2.108759)
  
  fig_FG_TYSR <- ggpredict(FG_mods_Bio[[i]][[lowest_AIC_model$model]], terms = c("logNumSp [1.3:3.6 by =0.2]", "Treatment", "Year1 [1:16, by = 5]")) %>% 
    data.frame() %>%
    mutate(facet = as.numeric(paste(facet)) + 2008) %>%
    ggplot(aes(y = predicted, x = exp(x), color = group, fill = group)) +
    theme(legend.position = "bottom")+
    facet_grid(.~facet) +
    scale_color_manual(breaks = mybreaks,
                       values = myvalues,
                       labels = mylabels) +
    scale_fill_manual(breaks = mybreaks,
                      values = myvalues,
                      labels = mylabels) +
    scale_x_continuous(breaks = c(1, 4, 16, 32)) +
    labs(title = "Biodiversity experiment",
         y = traitlabels_unit[[i]],
         x = "Species richness",
         color = "Treatment",
         fill = "Treatment") +
    geom_line() #+
  # geom_label(data = cbind.data.frame(label = (anova_Bio %>% filter(trait %in% i & Var %in% "Treatment:logNumSp:Year1"))$print,
  #                                    x = -Inf,
  #                                    y = Inf,
  #                                    facet = 2009),
  #            aes(y = y, x = x, label = label),
  #            color = "black",
  #            # x = -Inf,
  #            hjust = 0,
  #            # y = Inf,
  #            vjust = 1,
  #            label.size = 0,
  #            fill = NA,
  #            size = 2.108759)
  
  plots_FG_Bio[[i]] <- list(Treatment_Year = fig_FG_TY,
                            SR_Year = fig_FG_SRY,
                            Treatment_SR = fig_FG_SRT,
                            Treatment_SR_Year = fig_FG_TYSR)

}
  
rm(df, lin_mod_b, lin_mod_b_L, lin_mod_b_F, lin_mod_b_C3, lin_mod_b_C4, t1, t2,
   i, fig_TY, fig_SRY, fig_SRT, fig_TYSR, fig_FG_TY, fig_FG_SRY, fig_FG_SRT, 
   fig_FG_TYSR, AIC_header, lowest_AIC_model)


#### -> TABLE S1 ####
all_anova_for_print_Bio <- anova_Bio %>%
  mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
  mutate(P_print = gsub(P_print,
                        pattern = '<span style=" font-weight: bold; " >0</span>',
                        replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
  select(Var_print, trait, X2_print, P_print) %>%
  pivot_wider(id_cols = Var_print,
              names_from = trait,
              values_from = c(X2_print, P_print)) %>%
  select(Var_print,
         X2_print_Nmass,        P_print_Nmass,
         X2_print_DiasporeMass, P_print_DiasporeMass,
         X2_print_PlantHeight,  P_print_PlantHeight,
         X2_print_LeafArea,     P_print_LeafArea,
         X2_print_LMA,          P_print_LMA,
         X2_print_LDMC,         P_print_LDMC) %>%
  kbl(escape = F, caption  = "Biodiversity Experiment", col.names = c("Term", rep(c("X2", "P-value"), 6)))  %>%
  add_header_above(c(" ", "N content" = 2, "Seed Mass" = 2, "Plant Height" = 2, "Leaf Area" = 2, "LMA" = 2, "LDMC" = 2)) %>%
  kable_minimal(full_width = F, html_font = "Cambria") %>%
  kable_paper()


#___________________________________________________________________________####

## MODELS Oldfield Experiment ####
mods_OF <- list()
FG_mods_OF <- list()
anova_OF <- data.frame(Var = c(),
                        trait = c(), 
                        X2    = c(),
                        DF    = c(),
                        p     = c(),
                        P_print   = c(),
                        X2_print  = c(),
                        stars     = c(),
                        Var_print = c())

anova_for_print_OF <- list()
anova_for_print_FG_OF <- list()
mod_summary_OF <- list()
treatment_effects_OF <- list()
plots_OF <- list()
plots_FG_OF <- list()

for (i in c("Nmass", "DiasporeMass", "PlantHeight", "LeafArea", "LMA", "LDMC", "invsimpson", "PC1", "PC2")){
  
  ### data ####
  df <- CWM_OF %>% 
    merge(.,
          plot_data1 %>% 
            filter(experiment %in% "Oldfield") %>% 
            select(Plot, Treatment, Year, PC1, PC2) %>%
            rename("CWM_PC1_Biom"= "PC1",
                   "CWM_PC2_Biom" = "PC2"),
          by = c("Plot", "Treatment", "Year"), all.x = T) %>% 
    select(Block, Plot, Treatment, Year, Year1, C3_biom, C4_biom, F_biom, L_biom, contains(i)) %>%
    rename("y_b" = paste("CWM", i, "Biom", sep = "_")) %>%
     mutate(Treatment=`contrasts<-`(factor(Treatment), , OFc1))
  
  
  ### models ####
  lin_mod_b = lmer(y_b ~ Treatment * Year1                + (1|Block/Plot) + (1|Year), data = df)
  mods_OF[[i]] <- lin_mod_b
  
  lin_mod_b_L  <- update(lin_mod_b, .~. + L_biom)
  lin_mod_b_F  <- update(lin_mod_b, .~. + F_biom)
  lin_mod_b_C3 <- update(lin_mod_b, .~. + C3_biom)
  lin_mod_b_C4 <- update(lin_mod_b, .~. + C4_biom)
  
  FG_mods_OF[[i]] <- list(base_model   = lin_mod_b,
                           legume_model = lin_mod_b_L,
                           forb_model   = lin_mod_b_F,
                           C3_model     = lin_mod_b_C3,
                           C4_model     = lin_mod_b_C4)
  
  ### anova table (df & for print) ####
  anova_OF <- rbind.data.frame(anova_OF,
                                car::Anova(lin_mod_b) %>% 
                                  data.frame() %>%
                                  rename("X2" = "Chisq",
                                         "DF" = "Df",
                                         "p" = "Pr..Chisq.") %>%
                                  rownames_to_column("Var") %>%
                                  mutate(trait = i)  %>%
                                  mutate(P_print2 = case_when(p < 0.001 ~ "p < 0.001", .default = paste("p =", round(p, digits = 3))),
                                         P_print = case_when(p < 0.001 ~ "<0.001", .default = paste(round(p, digits = 3))),
                                         X2_print = round(X2, digits = 4),
                                         stars = case_when(p < 0.001 ~ "***",
                                                           p > 0.001 & p <0.01 ~ "**",
                                                           p > 0.01  & p < 0.05 ~ "*",
                                                           p > 0.05 & p <0.1 ~ ".",
                                                           .default = ""),
                                         Var_print = case_when(Var %in% "Treatment" ~ "Treatment (T)",
                                                               Var %in% "logNumSp"  ~ "log Species Richness (SR)",
                                                               Var %in% "Year1"     ~ "Year (Y)",
                                                               Var %in% "I(Year1^2)" ~ "Year^2 (Y2)",
                                                               Var %in% "Treatment:logNumSp"   ~ "T x SR",
                                                               Var %in% "Treatment:Year1"      ~ "T x Y",
                                                               Var %in% "Treatment:I(Year1^2)" ~ "T x Y2",
                                                               Var %in% "logNumSp:Year1"       ~ "SR x Y",
                                                               Var %in% "logNumSp:I(Year1^2)"  ~ "SR x Y2",
                                                               Var %in% "Treatment:logNumSp:Year1"      ~ "T x SR x Y",
                                                               Var %in% "Treatment:logNumSp:I(Year1^2)" ~ "T x SR x Y2"),
                                         print = paste(Var_print, ":", P_print2, stars)) 
  ) 
  
  
  anova_for_print_OF[[i]] <- anova_OF %>%
    filter(trait %in% i) %>%
    select(!trait) %>%
    mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
    mutate(P_print = gsub(P_print,
                          pattern = '<span style=" font-weight: bold; " >0</span>',
                          replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
    select(Var_print, X2_print, DF, P_print) %>%
    kbl(escape = F, 
        caption = paste("Biodiversity experiment:", traitlabels[[i]], sep = " ") , 
        col.names = c("Term", "X2", "DF",  "P-value"))  %>%
    kable_minimal(full_width = F, html_font = "Cambria") %>%
    kable_paper()
  
  
  
  
  #### -> Tables S4,6,8,10,12,14,16 ####
  AIC_header <-  rbind.data.frame("AIC", 
                                  AIC(lin_mod_b, lin_mod_b_L, lin_mod_b_F, lin_mod_b_C3, lin_mod_b_C4) %>%
                                    round(0) %>%
                                    data.frame %>%
                                    select(AIC)) %>%
    cbind.data.frame(c(1, rep(2, 5)))
  
  anova_for_print_FG_OF[[i]] <- car::Anova(lin_mod_b) %>% 
    data.frame() %>%
    rownames_to_column("Var") %>%
    mutate(model = "base_model") %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_L) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "legume_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_F) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "forb_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_C3) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "C3_model")) %>%
    
    rbind.data.frame(
      car::Anova(lin_mod_b_C4) %>% 
        data.frame() %>%
        rownames_to_column("Var") %>%
        mutate(model = "C4_model")) %>%
    
    rename("X2" = "Chisq",
           "p" = "Pr..Chisq.") %>%
    mutate(Var = case_when(Var %in% c("L_biom", "F_biom", "C3_biom", "C4_biom") ~  "FG_biomass", .default = Var)) %>%
    mutate(P_print = case_when(p < 0.001 ~ "<0.001", .default = paste(round(p, digits = 3))),
           X2_print = round(X2, digits = 4),
           Var_print = case_when(Var %in% "Treatment" ~ "Treatment (T)",
                                 Var %in% "Year1"     ~ "Year (Y)",
                                 Var %in% "Treatment:Year1"      ~ "T x Y",
                                 Var %in% "FG_biomass" ~ "FG biomass")) %>%
    mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
    mutate(P_print = gsub(P_print,
                          pattern = '<span style=" font-weight: bold; " >0</span>',
                          replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
    select(Var_print, model, X2_print, P_print) %>%
    
    pivot_wider(id_cols = Var_print,
                values_from = c(X2_print, P_print),
                names_from = model) %>%
    select(Var_print, 
           X2_print_base_model,  P_print_base_model,
           X2_print_legume_model, P_print_legume_model,              
           X2_print_forb_model,  P_print_forb_model,
           X2_print_C3_model, P_print_C3_model,
           X2_print_C4_model, P_print_C4_model)  %>%
    
    kbl(escape = F, 
        caption = paste("Oldfield experiment:", traitlabels[[i]], sep = " "),
        col.names = c("Term", rep(c("X2", "P-value"), 5)))  %>%
    add_header_above(AIC_header) %>%
    add_header_above(c(" ", "Base model" = 2, "Legume model" = 2, "Forb model" = 2, "C3 model" = 2, "C4 model" = 2))%>%
    kable_paper()
  
  ### model summary for print ####
  t1 <- lin_mod_b %>% 
    summary() %>%
    coef() %>% 
    data.frame() %>%
    round(digits = 4) %>%
    rename("p" = "Pr...t..") %>%
    mutate(p = cell_spec(p, bold = ifelse(p < 0.05, TRUE, FALSE))) %>%
    mutate(p = gsub(p, pattern = '<span style=" font-weight: bold; " >0</span>', replacement = '<span style=" font-weight: bold;    " ><0.001</span>')
    ) %>%
    rename("SE" = "Std..Error",
           "DF" = "df",
           "t-value" = "t.value",
           "p-value" = "p") %>%
    kbl(escape = F, 
        caption = paste("Biodiversity experiment: CWM", i, sep = " ")) %>%
    kable_minimal(full_width = F, html_font = "Cambria") %>%
    kable_paper()
  
  t2 <- VarCorr(lin_mod_b) %>%
    as.data.frame() %>%
    select(grp, vcov, sdcor) %>%
    mutate(vcov = round(vcov, digits=4),
           sdcor = round(sdcor, digits = 4)) %>%
    rename("Random effect" = "grp",
           "Variance" = "vcov",
           "DF" = "sdcor") %>%
    kbl(escape = F) %>%
    kable_minimal(full_width = F,  html_font = "Cambria") %>%
    kable_paper()
  
  mod_summary_OF[[i]] <- list(coefs = t1,
                              ranefs = t2)
  
  
  ### emmeans pairwise ####
  # treatment_effects_OF[[i]] <- (emmeans(lin_mod_b,  pairwise ~ Treatment))$contrasts %>%
  #   data.frame()
  
  ### plots ####
  #### -> Figure S5-7,10a ####
  fig_TY <- expand_grid(Year1 = c(1:16),
                        Treatment = c("Control",
                                      "Fenced",
                                      "Insecticide",
                                      "SoilDrenchFungicide",
                                      "FoliarFungicide",
                                      "AllPesticides")) %>%
    mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
           Year = Year1+2008) %>%
    mutate(fixed_lb = predict(mods_OF[[i]],  newdata = . , re.form = NA)) %>%
    
    merge(
      expand_grid(Year1 = c(1:16),
                  Treatment = c("Control",
                                "Fenced",
                                "Insecticide",
                                "SoilDrenchFungicide",
                                "FoliarFungicide",
                                "AllPesticides")) %>%
        mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
               Year = Year1+2008) %>%
        filter(!Year %in% c(2019, 2020)) %>%
        mutate(full_lb = predict(mods_OF[[i]],  newdata = . , re.form = ~(1|Year))),
      
      all = T, by = c("Year1", "Treatment", "logNumSp", "Year")
    ) %>%
    
    ggplot(aes(y = fixed_lb, x = Year, color = Treatment)) +
    labs(title = "Oldfield experiment",
         y = traitlabels_unit[[i]] )+
    scale_color_manual(breaks = mybreaks, labels= mylabels, values= myvalues)+
    scale_x_continuous(breaks = ylabels) +
    geom_line() +
    geom_line(aes(y=full_lb), linetype="dotted") +
    geom_label(label = (anova_OF %>%
                          filter(trait %in% i) %>%
                          filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
                          select(print) %>%
                          summarize(print = paste(print, collapse = "\n")))$print,
               color = "black",
               x = -Inf, hjust = 0,
               y = Inf, vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)
  
  
  plots_OF[[i]] <- fig_TY
  
  
  ### plots of functional group model with lowest AIC ####
  
  #### -> Figure S5-7,10b ####
  
  # find functional group model with lowest AIC
  lowest_AIC_model <- cbind.data.frame(
    model = c("legume_model", "forb_model", "C3_model", "C4_model"),
    var = c("L_biom", "F_biom", "C3_biom", "C4_biom"),
    AIC =  (AIC(FG_mods_OF[[i]]$legume_model,
                FG_mods_OF[[i]]$forb_model,
                FG_mods_OF[[i]]$C3_model,
                FG_mods_OF[[i]]$C4_model))$AIC) %>%
    filter(AIC %in% min(AIC))
  
  
  ### plots
  fig_FG_TY <- expand_grid(Year1 = c(1:16),
                         Treatment = c("Control",
                                       "Fenced",
                                       "Insecticide",
                                       "SoilDrenchFungicide",
                                       "FoliarFungicide",
                                       "AllPesticides"),
                        !!lowest_AIC_model$var := mean(CWM_OF[,lowest_AIC_model$var], na.rm =T)
  ) %>%
    mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
           Year = Year1+2008) %>%
    mutate(fixed_lb = predict(FG_mods_OF[[i]][[lowest_AIC_model$model]],  newdata = . , re.form = NA)) %>%
    
    merge(
      expand_grid(Year1 = c(1:16),
                  Treatment = c("Control",
                                "Fenced",
                                "Insecticide",
                                "SoilDrenchFungicide",
                                "FoliarFungicide",
                                "AllPesticides"),
                  !!lowest_AIC_model$var := mean(CWM_OF[,lowest_AIC_model$var], na.rm =T)) %>%
        mutate(logNumSp = mean(log(c(1, 4, 16, 32))),
               Year = Year1+2008) %>%
        filter(!Year %in% c(2019, 2020)) %>%
        mutate(full_lb = predict(FG_mods_OF[[i]][[lowest_AIC_model$model]],  newdata = . , re.form = ~(1|Year))),
      
      all = T, by = c("Year1", "Treatment", "logNumSp", "Year")
    ) %>%
    
    ggplot(aes(y = fixed_lb, x = Year, color = Treatment)) +
    labs(title = paste("Oldfield experiment,", gsub(lowest_AIC_model$model, pattern = "_", replacement = " ")),
         y = traitlabels_unit[[i]] )+
    scale_color_manual(breaks = mybreaks, labels= mylabels, values= myvalues)+
    scale_x_continuous(breaks = ylabels) +
    geom_line() +
    geom_line(aes(y=full_lb), linetype="dotted") 
    # geom_label(label = (anova_OF %>%
    #                       filter(trait %in% i) %>%
    #                       filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
    #                       select(print) %>%
    #                       summarize(print = paste(print, collapse = "\n")))$print,
    #            color = "black",
    #            x = -Inf, hjust = 0,
    #            y = Inf, vjust = 1,
    #            label.size = 0,
    #            fill = NA,
    #            size = 2.108759)
  
  
  plots_FG_OF[[i]] <- fig_FG_TY
  
}

rm(df, lin_mod_b, lin_mod_b_L, lin_mod_b_F, lin_mod_b_C3, lin_mod_b_C4, t1, t2,fig_FG_TY,
   i, fig_TY, AIC_header, lowest_AIC_model)

#### -> Table S2 ####
all_anova_for_print_OF <- anova_OF %>%
  mutate(P_print = cell_spec(P_print, bold = ifelse(p < 0.05, TRUE, FALSE)) ) %>%
  mutate(P_print = gsub(P_print,
                        pattern = '<span style=" font-weight: bold; " >0</span>',
                        replacement = '<span style=" font-weight: bold;    " ><0.001</span>')) %>%
  select(Var_print, trait, X2_print, P_print) %>%
  pivot_wider(id_cols = Var_print,
              names_from = trait,
              values_from = c(X2_print, P_print)) %>%
  select(Var_print,
         X2_print_Nmass,        P_print_Nmass,
         X2_print_DiasporeMass, P_print_DiasporeMass,
         X2_print_PlantHeight,  P_print_PlantHeight,
         X2_print_LeafArea,     P_print_LeafArea,
         X2_print_LMA,          P_print_LMA,
         X2_print_LDMC,         P_print_LDMC) %>%
  kbl(escape = F, caption  = "Oldfield Experiment", col.names = c("Term", rep(c("X2", "P-value"), 6)))  %>%
  add_header_above(c(" ", "N content" = 2, "Seed Mass" = 2, "Plant Height" = 2, "Leaf Area" = 2, "LMA" = 2, "LDMC" = 2)) %>%
  kable_minimal(full_width = F, html_font = "Cambria") %>%
  kable_paper()



#___________________________________________________________________________####

# POSTHOC TESTS & PLOTS ####

### emmeans ####
posthoc_Bio <- list()
for (i in c("Nmass", "LeafArea", "LMA", "LDMC", "PlantHeight", "DiasporeMass", "invsimpson", "PC1", "PC2")){
posthoc_Bio[[i]] <- (emmeans(mods_Bio[[i]], pairwise ~ Treatment|Year1, at = list(Year1 = seq(1,16, 0.2))))$contrasts
}

posthoc_OF <- list()
for (i in c("Nmass", "LeafArea", "LMA", "LDMC", "PlantHeight", "DiasporeMass", "invsimpson", "PC1", "PC2")){
  posthoc_OF[[i]] <- (emmeans(mods_OF[[i]], pairwise ~ Treatment|Year1, at = list(Year1 = seq(1,16, 0.2))))$contrasts
}


### Treatment effect ~ Time ####
plots_TrtEff_Bio <- list()
plots_TrtEff_OF <- list()
plots_TrtEff <- list()


for (i in c("Nmass", "DiasporeMass", "PlantHeight", "LeafArea", "LMA", "LDMC", "invsimpson", "PC1", "PC2")){
  df_Bio <- posthoc_Bio[i] %>% 
    data.frame() %>% 
    rename_with(~str_remove(., paste(i, ".", sep = "")))%>%
    filter(grepl(x = contrast, pattern = "Control")) %>%
    mutate(Year = Year1 + 2008,
           Treatment_effect = estimate * (-1),
           Significance = case_when(p.value < 0.05 ~ "p < 0.05", .default ="p > 0.05"),
           Treatment = gsub(contrast, pattern = "Control - ", replacement = ""))
  
  df_OF <- posthoc_OF[i] %>% 
    data.frame() %>% 
    rename_with(~str_remove(., paste(i, ".", sep = "")))%>%
    filter(grepl(x = contrast, pattern = "Control")) %>%
    mutate(Year = Year1 + 2008,
           Treatment_effect = estimate * (-1),
           Significance = case_when(p.value < 0.05 ~ "p < 0.05", .default ="p > 0.05"),
           Treatment = gsub(contrast, pattern = "Control - ", replacement = "")) 
  
  
  #### -> Figure S8 ####
  plots_TrtEff[[i]] <- df_Bio %>%
    mutate(experiment = "Biodiversity") %>%
    rbind.data.frame(df_OF %>% 
                       mutate(experiment = "Oldfield")) %>%
    
    ggplot(aes(y = Treatment_effect, x = Year)) +
    theme(axis.title.x = element_blank())+
    facet_grid(.~experiment)+
    scale_color_manual(breaks = mybreaks, labels = mylabels, values=myvalues) +
    scale_x_continuous(breaks = ylabels) +
    scale_linetype_manual(breaks = sigbreaks, values = sigvalues) +
    labs(y = traitlabels_unit_2[[i]]) +
    geom_line(aes(linetype = Significance, color = Treatment)) +
    geom_abline(intercept = 0, slope = 0, color = "grey") +
    geom_abline(intercept = 0, slope = 0, color = "grey") +
    geom_label(data =   anova_Bio %>%
                 filter(trait %in% i) %>%
                 filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
                 select(print) %>%
                 summarize(print = paste(print, collapse = "\n")) %>%
                 mutate(experiment = "Biodiversity") %>% 
                 
                 rbind.data.frame(
                   anova_OF %>%
                     filter(trait %in% i) %>%
                     filter(Var %in% c("Treatment", "Year1", "Treatment:Year1")) %>%
                     select(print) %>%
                     summarize(print = paste(print, collapse = "\n")) %>%
                     mutate(experiment = "Oldfield")    ),
               
               aes(label = print),
               color = "black",
               # x = -Inf, hjust = 0,
               x = 2012, hjust = 0,
               y = Inf, vjust = 1,
               label.size = 0,
               fill = NA,
               size = 2.108759)

  }

rm(i, df_Bio, df_OF)

### Correlation between early and late treatment effects ####
rect_helper <- list(Nmass = Inf,
                    LeafArea = Inf,
                    LMA = -Inf,
                    LDMC = -Inf,
                    PlantHeight = Inf,
                    DiasporeMass = Inf,
                    invsimpson = -Inf)

plots_effect_correlations <- list()
plots_effect_correlations2 <- list()
for (i in c("Nmass", "LeafArea", "LMA", "LDMC", "PlantHeight", "DiasporeMass", "invsimpson", "PC1", "PC2")){
  df <- posthoc_Bio[i] %>% 
    data.frame() %>% 
    rename_with(~str_replace(., paste(i, ".", sep = ""), "Biodiversity."))%>%
    rename("contrast" = "Biodiversity.contrast",
           "Year1" = "Biodiversity.Year1") %>%
    filter(grepl(x = contrast, pattern = "Control")) %>%
    filter(Year1 %in% c(1,16)) %>%
    mutate(Biodiversity.Treatment_effect = Biodiversity.estimate * (-1),
           Treatment = gsub(contrast, pattern = "Control - ", replacement = ""),
           Biodiversity.Significance = case_when(Biodiversity.p.value < 0.05 ~ "p < 0.05", .default ="p > 0.05")) %>%
    
    merge(
      posthoc_OF[i] %>% 
        data.frame() %>% 
        rename_with(~str_replace(., paste(i, ".", sep = ""), "Oldfield."))%>%
        rename("contrast" = "Oldfield.contrast",
               "Year1" = "Oldfield.Year1") %>%
        filter(grepl(x = contrast, pattern = "Control")) %>%
        filter(Year1 %in% c(1, 16)) %>%
        mutate(Oldfield.Treatment_effect = Oldfield.estimate * (-1),
               Treatment = gsub(contrast, pattern = "Control - ", replacement = ""),
               Oldfield.Significance = case_when(Oldfield.p.value < 0.05 ~ "p < 0.05", .default ="p > 0.05")),
      
      by = c("contrast", "Year1", "Treatment")
    ) %>%
    mutate(Year = Year1+2008)
  
  #### -> Figure 2 ####
  plots_effect_correlations[[i]] <- df %>%
    filter(Year1 %in% 16) %>%
    
    ggplot(aes(y = Oldfield.Treatment_effect, 
               x = Biodiversity.Treatment_effect)) +
      labs(y = "Treatment effect \u00B1 SE\nin the Oldfield experiment",
           x = "Treatment effect \u00B1 SE\nin the Biodiversity experiment",
           title = traitlabels[[i]]) +
      geom_rect(ymin = 0, xmin = 0, ymax =  (-1)*(rect_helper[[i]]), xmax = rect_helper[[i]], fill = "#FAFAFA") +
      geom_rect(ymin = 0, xmin = 0, ymax =  rect_helper[[i]], xmax = (-1)*(rect_helper[[i]]), fill = "#FAFAFA") +
      geom_rect(ymin = 0, xmin = 0, ymax =  rect_helper[[i]], xmax = rect_helper[[i]], fill = "#F5F5F5") +
      geom_hline(yintercept = 0, color = "grey", linetype = "dotted")+
      geom_vline(xintercept = 0, color = "grey", linetype = "dotted") +
      # geom_abline(intercept = 0, slope = 1, color = "grey") +
      # geom_smooth(method = "lm", color = "grey", alpha = 0.2, se = F)+
      scale_color_manual(breaks = mybreaks, values=myvalues, labels = mylabels) +
      scale_linetype_manual(breaks = sigbreaks, values = sigvalues) +
      scale_size_manual(breaks = sigbreaks, values = sigvalues2) +
      geom_point(aes(color = Treatment)) +
      geom_errorbar(aes(ymin = Oldfield.Treatment_effect-Oldfield.SE, 
                        ymax = Oldfield.Treatment_effect+Oldfield.SE,
                        color = Treatment,
                        linetype = Oldfield.Significance,
                        size = Oldfield.Significance)) +
      geom_errorbarh(aes(xmin = Biodiversity.Treatment_effect- Biodiversity.SE,
                         xmax = Biodiversity.Treatment_effect+ Biodiversity.SE,
                         color = Treatment,
                         linetype = Biodiversity.Significance,
                         size = Biodiversity.Significance)) +
      stat_cor(label.x.npc = "left", label.y.npc = "bottom")
  
  
}

rm(i, df, rect_helper)





##__________________________________________________________________________####

# SPECIES CHANGE ####

### heatmap ####
shared_species <- sp_change %>%
  filter(exp %in% c("Biodiversity experiment", "Oldfield experiment")) %>%
  group_by(Species, Functional.group) %>%
  summarize(n_exp = n_distinct(exp), .groups = "drop") %>%
  filter(n_exp == 2) 

# species_faces <- setNames(
#   ifelse(unique(sp_change$Species) %in% shared_species, "bold.italic", "italic"),
#   unique(sp_change$Species)
# )

#### -> Figure 3 ####
heatmap <-
  sp_change %>%
  mutate(stars = case_when(p.value > 0.1 ~ "",
                           p.value > 0.5 ~ ".",
                           p.value > 0.01 ~ "*",
                           p.value > 0.001 ~ "**",
                           p.value < 0.001 ~ "***",
                           .default = "fuck"),
         Species = case_when(Species %in% "Achillea millefolium (lanulosa)" ~ "Achillea millefolium",
                             Species %in%"Ambrosia artemisiifolia elatior" ~ "Ambrosia artemisiifolia", 
                             .default=Species)) %>%
  arrange(Functional.group, Species) %>%
  mutate(Species = factor(Species, unique(Species)),
         Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides"))) %>% 
  # filter(p.value > 0.05) %>%
  filter(!Species %in% c("Rhus glabra", "Rosa arkansana")) %>%
  mutate(Species = ifelse(Species %in% shared_species$Species, glue("<b><i>{Species}</i></b>"), glue("<i>{Species}</i>"))) %>%
  
  ggplot(aes(y = Species, x = Treatment, fill = Year2.trend)) +
  facet_grid(Functional.group~exp, scales = "free", space = "free",
             labeller = labeller(Functional.group = c(C3 = "C3\nGrasses", C4 = "C4\nGrasses", 'F' = "Non-leguminous forbs", L = "Legumes" ))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = ggtext::element_markdown(),
        # axis.text.y = element_text(face = "italic"),
        axis.title = element_blank(),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        panel.border = element_blank(),
        legend.key.width = unit(0.6, "cm"),
        legend.position = "bottom",
        legend.justification = "center") +
  scale_x_discrete(breaks = c("Control", "Fenced", "Insecticide", "FoliarFungicide", 
                              "SoilDrenchFungicide", "AllPesticides"),
                   labels = c("Control", "Fenced", "Insecticide", "Foliar fungicide", 
                              "Soil drench fungicide", "All pesticides")) +
  scale_y_discrete(limits=rev) +
  labs(fill = "Change over time")+
  geom_tile() +
  # geom_hline(yintercept = c(3.5, 8.5, 18.5)) +
  scale_fill_gradient2(low = '#D41159', mid = 'lightgrey', high = '#1A85FF')+
  geom_text(aes(label = stars))

### Treatment effect ####
#### -> Figure S12 ####
plot_sp_TrtEff_cor <-
  sp_TrtEff %>%
  filter(!grepl("Fenced", contrast)) %>%
  mutate(Treatment = gsub(contrast, pattern = "Control - ", replacement = "") %>%
           factor(levels = c("Insecticide",
                             "SoilDrenchFungicide",
                             "FoliarFungicide",
                             "AllPesticides")),
         Significance = case_when(p.value < 0.05 ~ "p < 0.05", .default = "p > 0.05"),
         exp = gsub(exp, pattern = " experiment", replacement = "")) %>%
  filter(!Functional.group %in% "W") %>%
  filter(!contrast %in% NA) %>%
  pivot_wider(id_cols = c(Species, Functional.group, Treatment),
              values_from = c(estimate, SE, Significance),
              names_from = exp) %>%
  filter(!(estimate_Oldfield %in% NA)) %>%
  filter(!(estimate_Biodiversity %in% NA)) %>%
  
  ggplot(aes(y = estimate_Oldfield, x = estimate_Biodiversity, color = Treatment)) +
  facet_wrap(~Treatment,
             labeller = labeller(Treatment = c("Insecticide" = "Insecticide",
                                               "SoilDrenchFungicide" = "Soil drench fungicide",
                                               "FoliarFungicide" = "Foliar fungicide",
                                               "AllPesticides" = "All pesticides"))) +
  labs(shape = "Functional group",
       y = "Treatment effect on species biomass\nin the Oldfield experiment",
       x = "Treatment effect on species biomass\nin the Biodiversity experiment") +
  guides(color = FALSE)+
  scale_linetype_manual(breaks = sigbreaks, values = sigvalues) +
  scale_color_manual(breaks = mybreaks, labels = mylabels, values = myvalues) +
  scale_shape_manual(breaks = c("C3", "C4", "F", "L"),
                     labels = c("C3 grasses", "C4 grasses", "Non-leguminous forbs", "Legumes"),
                     values = c(15, 16, 17, 8))+
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") + 
  geom_text_repel(aes(label = Species),
                  size = 3,
                  segment.color = 'grey',
                  segment.size = 0.2,
                  color = "grey50",
                  force = 80
                  ) +
  geom_point(aes(shape = Functional.group)) +
  # geom_errorbar(aes(ymin = estimate_Oldfield - SE_Oldfield,
  #                   ymax = estimate_Oldfield + SE_Oldfield,
  #                   linetype = Significance_Biodiversity)) +
  # geom_errorbarh(aes(xmin = estimate_Biodiversity - SE_Biodiversity,
  #                    xmax = estimate_Biodiversity + SE_Biodiversity,
  #                    linetype = Significance_Biodiversity)) +
  stat_cor(color = "black", aes(group = "Treatment")) 




df <- sp_TrtEff %>%
    filter(!contrast %in% NA) %>%
    mutate(Treatment = gsub(contrast, pattern = "Control - ", replacement = "") %>%
             factor(levels = c("Fenced", 
                               "Insecticide",
                               "SoilDrenchFungicide",
                               "FoliarFungicide",
                               "AllPesticides")),
           Significance = case_when(p.value < 0.05 ~ "p < 0.05", .default = "p > 0.05")) %>%
  pivot_longer(cols = c(Nmass, LeafArea, LMA, PlantHeight, DiasporeMass, LDMC),
               names_to = "trait",
               values_to = "trait_value") %>%
  mutate(trait = factor(trait, levels = c("Nmass", "DiasporeMass", "PlantHeight",
                                          "LeafArea", "LMA", "LDMC"))) %>%
  filter(!is.na(trait_value))

df_sp_labs <- df %>%
  filter(Species %in% c("Monarda fistulosa",
                       "Poa pratensis",
                       "Liatis aspera",
                       "Solidago nemoralis",
                       "Lespedeza capitata",
                       "Koeleria cristata",
                       "Sorghastrum nutans",
                       (sp_TrtEff %>% slice(which.min(Nmass)))$Species,
                       (sp_TrtEff %>% slice(which.max(Nmass)))$Species,
                       (sp_TrtEff %>% slice(which.min(DiasporeMass)))$Species,
                       (sp_TrtEff %>% slice(which.max(DiasporeMass)))$Species,
                       (sp_TrtEff %>% slice(which.min(PlantHeight)))$Species,
                       (sp_TrtEff %>% slice(which.max(PlantHeight)))$Species,
                       (sp_TrtEff %>% slice(which.min(LeafArea)))$Species,
                       (sp_TrtEff %>% slice(which.max(LeafArea)))$Species,
                       (sp_TrtEff %>% slice(which.min(LMA)))$Species,
                       (sp_TrtEff %>% slice(which.max(LMA)))$Species,
                       (sp_TrtEff %>% slice(which.min(LDMC)))$Species,
                       (sp_TrtEff %>% slice(which.max(LDMC)))$Species)) %>%
  mutate(y_min = min(estimate, na.rm=T),
         y_offset = diff(range(estimate)),
         y_below_axis = y_min - 0.6*y_offset) %>%
  select(Species, trait, trait_value, y_min, y_offset, y_below_axis) %>%
  unique() %>%
  filter(!is.na(trait_value)) %>%
  mutate(exp = "Oldfield experiment",
         species_label = paste("italic(\"", Species, "\")", sep = "")) %>%
  filter(Species %in% c("Monarda fistulosa",
                        "Poa pratensis",
                        "Liatis aspera",
                        "Solidago nemoralis",
                        "Lespedeza capitata",
                        "Koeleria cristata",
                        "Sorghastrum nutans") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(Nmass)))$Species,
                           (sp_TrtEff %>% slice(which.max(Nmass)))$Species) &
              trait %in% "Nmass") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(DiasporeMass)))$Species,
                           (sp_TrtEff %>% slice(which.max(DiasporeMass)))$Species) &
              trait %in% "DiasporeMass") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(PlantHeight)))$Species,
                           (sp_TrtEff %>% slice(which.max(PlantHeight)))$Species) &
              trait %in% "PlantHeight") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(LeafArea)))$Species,
                           (sp_TrtEff %>% slice(which.max(LeafArea)))$Species) &
              trait %in% "LeafArea") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(LMA)))$Species,
                           (sp_TrtEff %>% slice(which.max(LMA)))$Species) &
              trait %in% "LMA") |
           (Species %in% c((sp_TrtEff %>% slice(which.min(LDMC)))$Species,
                           (sp_TrtEff %>% slice(which.max(LDMC)))$Species) &
              trait %in% "LDMC")  )



#### -> Figure S11 ####
plot_sp_TrtEff2 <-
    
    ggplot(df, aes(y = estimate, x = trait_value, color = Treatment)) +
      theme(legend.position = "bottom") +
    coord_cartesian(clip = "off") +
    facet_grid(exp~trait, scales = "free_x",
               labeller=labeller(trait = c("Nmass" = "Tissue N",
                                           "DiasporeMass" = "Seed mass",
                                           "PlantHeight" = "Plant height",
                                           "LeafArea" = "Leaf area",
                                           "LMA" = "LMA",
                                           "LDMC" = "LDMC"))) +
    labs(x = "Trait value",
         shape = "Functional group",
         y = "Treatment effect \u00B1 SE") +
    guides(linetype = FALSE, color = guide_legend(override.aes = list(size = 2)))+
    coord_cartesian(ylim = c(-1, 1))+
    scale_color_manual(breaks = mybreaks, labels = mylabels, values = myvalues) +
    scale_shape_manual(breaks = c("C3", "C4", "F", "L"),
                       labels = c("C3 grasses", "C4 grasses", "Non-leguminous forbs", "Legumes"),
                       values = c(15, 16, 17, 8))+
    geom_hline(yintercept = 0, color = "grey")  + 
    geom_point(aes(shape = Functional.group), size=2) +
    geom_smooth(method = "lm", se = F) +
    stat_cor() +
    geom_text_repel(data = df_sp_labs,
            aes(x = trait_value, y = y_below_axis, label = species_label),
            angle = 90,
            hjust = 0,
            vjust = 0.5,
            size = 3,
            color = "black",
            direction = "x",
            segment.color = "grey",
            parse=T)

rm(df, df_sp_labs)


    
#___________________________________________________________________________####
## DATA QUALITY ####

### get species level data ####
plotinfo_BigBio <- read.csv(paste(EnemyRemovalDIR, "data-raw/BigBio_All_Plots_info.csv", sep = "")) %>%
  select(!c(Exp)) 

plotinfo_Bio <- read.csv(paste(EnemyRemovalDIR, "/data-raw/E244Data_PlanFilev1.csv", sep = ""))

design <- plotinfo_BigBio %>%
  pivot_longer(cols = !c(Plot, NumSp, FgNum, Fgset, C3, C4, Forb, Legume, Woody),
               names_to = "Species",
               values_to = "Sown") %>%
  filter(Sown %in% 1) %>%
  mutate(Species_Diaz = case_when(Species %in% "Achillea.millefolium.lanulosa." ~ "Achillea millefolium",
                                  Species %in% "Agropyron.smithii"              ~ "Agropyron smithii",
                                  Species %in% "Amorpha.canescens"              ~ "Amorpha canescens",
                                  Species %in% "Andropogon.gerardi"             ~ "Andropogon gerardi",
                                  Species %in% "Asclepias.tuberosa"             ~ "Asclepias tuberosa",
                                  Species %in% "Aster.azureus"                  ~ "Symphyotrichum oolentangiense",
                                  Species %in% "Astragalus.canadensis"          ~ "Astragalus canadensis",
                                  Species %in% "Baptisia.leucantha"             ~ "Baptisia alba",
                                  Species %in% "Bouteloua.curtipendula"         ~ "Bouteloua curtipendula",
                                  Species %in% "Bouteloua.gracilis"             ~ "Bouteloua gracilis",
                                  Species %in% "Bromus.inermis"                 ~ "Bromus inermis",
                                  Species %in% "Buchloe.dactyloides"            ~ "Buchloe dactyloides",
                                  Species %in% "Calamagrostis.canadensis"       ~ "Calamagrostis.canadensis",
                                  Species %in% "Coreopsis.palmata"              ~ "Coreopsis palmata",
                                  Species %in% "Elymus.canadensis"              ~ "Elymus canadensis",
                                  Species %in% "Koeleria.cristata"              ~ "Koeleria cristata",
                                  Species %in% "Leersia.oryzoides"              ~ "Leersia oryzoides",
                                  Species %in% "Lespedeza.capitata"             ~ "Lespedeza capitata",
                                  Species %in% "Liatris.aspera"                 ~ "Liatris aspera",
                                  Species %in% "Lupinus.perennis"               ~ "Lupinus perennis",
                                  Species %in% "Monarda.fistulosa"              ~ "Monarda fistulosa",
                                  Species %in% "Panicum.virgatum"               ~ "Panicum virgatum",
                                  Species %in% "Petalostemum.candidum"          ~ "Petalostemum candidum", # no traits
                                  Species %in% "Petalostemum.purpureum"         ~ "Dalea purpurea",
                                  Species %in% "Petalostemum.villosum"          ~ "Dalea villosa",
                                  Species %in% "Poa.pratensis"                  ~ "Poa pratensis",
                                  Species %in% "Quercus.ellipsoidalis"          ~ "Quercus ellipsoidalis", # exclude
                                  Species %in% "Quercus.macrocarpa"             ~ "Quercus macrocarpa", # exclude
                                  Species %in% "Rudbeckia.serotina"             ~ "Echinacea serotina",
                                  Species %in% "Schizachyrium.scoparium"        ~ "Schizachyrium scoparium",
                                  Species %in% "Solidago.nemoralis"             ~ "Solidago nemoralis",
                                  Species %in% "Solidago.rigida"                ~ "Solidago rigida",
                                  Species %in% "Sorghastrum.nutans"             ~ "Sorghastrum nutans",
                                  Species %in% "Sporobolus.cryptandrus"         ~ "Sporobolus cryptandrus",
                                  Species %in% "Stipa.spartea"                  ~ "Stipa spartea",
                                  Species %in% "Trifolium.hybridum"             ~ "Trifolium hybridum",
                                  Species %in% "Vicia.villosa"                  ~ "Vicia villosa",
                                  Species %in% "Zizia.aurea"                    ~ "Zizia aurea"),
         Species2 = as.character(paste(Species)),
         Species2 = gsub(Species2, pattern = "\\.", replacement = " "),
         Species2 = gsub(Species2, pattern ="Achillea millefolium lanulosa ", replacement = "Achillea millefolium")) %>%
  merge(plotinfo_Bio %>% 
          select("BigBioPlot", "ER244.Plot"),
        by.x = "Plot", by.y = "BigBioPlot") %>%
  rename(BigBioPlot = Plot,
         Plot = ER244.Plot) %>%
  select(BigBioPlot, Plot, Species, Species2, Species_Diaz) %>%
  mutate(Plot_Species_Diaz = paste(Plot, Species_Diaz),
         Plot_Species = paste(Plot, Species2))

speciesinfo <- read.csv(paste(EnemyRemovalDIR, "/data-derived/speciesinfo_extended.csv", sep = ""))

biom_ER_Bio <- read.csv(paste(EnemyRemovalDIR, "/data-derived/E244_biomass_clean.csv", sep = "")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", 
                                       "Fenced",
                                       "Insecticide", 
                                       "SoilDrenchFungicide", 
                                       "FoliarFungicide",
                                       "AllPesticides"))) %>%
  merge(speciesinfo %>% select(Species, Family, Functional.group), by = "Species", all.x = T)

biom_ER_OF <- read.csv(paste(EnemyRemovalDIR, "/data-derived/E245_biomass_clean.csv", sep = "")) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Control", 
                                       "Fenced",
                                       "Insecticide", 
                                       "SoilDrenchFungicide", 
                                       "FoliarFungicide",
                                       "AllPesticides"))) %>%
  merge(speciesinfo %>% select(Species, Family, Functional.group), by = "Species", all.x = T)

DiazTraits <- read.csv(paste(EnemyRemovalDIR, "data-raw/CDR_traits_SDaiz.csv", sep = "")) 


# make species filters across both cover data
sp_species <- c("Achillea millefolium (lanulosa) ", 
                "Achillea millefolium (lanulosa)",
                "Agropyron repens",
                "Agropyron smithii",
                "Agrostis scabra",
                "Ambrosia artemisiifolia elatior",
                "Ambrosia coronopifolia",
                "Amorpha canescens",
                "Andropogon gerardi",
                "Anemone cylindrica",
                "Antennaria neglecta",
                "Arabis divaricarpa",
                "Aristida basiramea",
                "Aristida tuberculosa",
                "Artemisia (caudata) campestris",
                "Artemisia ludoviciana",
                "Asclepias syriaca",
                "Asclepias ovalifolia",
                "Asclepias tuberosa",
                "Aster azureus",
                "Astragalus canadensis",
                "Baptisia leucantha",
                "Berteroa incana",
                "Bouteloua curtipendula",
                "Bouteloua gracilis",
                "Bromus inermis",
                "Buchloe dactyloides",  
                "Calamagrostis canadensis",
                "Calamovilfa longifolia",
                "Campanula rotundifolia",
                "Chenopodium leptophyllum",
                "Coreopsis palmata",
                "Corylus americanus", 
                "Crepis tectorum",
                "Cyperus filiculmis",
                "Cyperus schweinitzii",
                "Danthonia spicata",
                "Delphinium virescens",
                "Desmodium canadense",
                "Digitaria ischaemum",
                "Elymus canadensis",
                "Equisetum arvense",
                "Equisetum hyemale",
                "Equisetum laevigatum",
                "Eragrostis spectabilis",
                "Erigeron canadensis",
                "Erigeron strigosus",
                "Eupatorium perfoliatum",
                "Euphorbia corollata",
                "Euphorbia glyptosperma",
                "Fragaria virginiana",
                "Hedeoma hispida",
                "Hieracium longipilum",
                "Koeleria cristata",
                "Lepidium densiflorum",
                "Leptoloma cognatum",
                "Lespedeza capitata",
                "Leersia oryzoides",        
                "Liatris aspera",
                "Linaria vulgaris",
                "Lupinus perennis",
                "Lychnis alba",
                "Mirabilis hirsuta",
                "Mollugo verticillata",
                "Monarda fistulosa",
                "Muhlenbergia racemosa",
                "Oenothera biennis",
                "Oxybaphus hirsutus",
                "Paspalum ciliatifolium",
                "Panicum oligosanthes",
                "Panicum perlongum",
                "Panicum virgatum",
                "Penstemon grandiflorus",          
                "Petalostemum candidum",
                "Petalostemum purpureum",
                "Petalostemum villosum",
                "Physalis heterophylla",
                "Physalis virginiana",
                "Poa pratensis",
                "Polygonatum canaliculatum",
                "Polygonum convolvulus",
                "Polygonum tenue",
                "Polygonella articulata",
                "Quercus alba",
                "Quercus ellipsoidalis",
                "Quercus macrocarpa",
                "Rosa arkansana",
                "Rubus alleghaniensis",
                "Rudbeckia serotina",
                "Rumex acetosella",
                "Schizachyrium scoparium",
                "Setaria lutescens (glauca)",
                "Silene antirrhina",
                "Smilacina stellata",
                "Solidago gigantea",
                "Solidago nemoralis",
                "Solidago rigida",
                "Sorghastrum nutans",
                "Sporobolus cryptandrus",
                "Stipa spartea",
                "Rhus typhina",
                "Scutellaria lateriflora",
                "Scutellaria parvula",
                "Smilacina racemosa",
                "Solanum americanum",
                "Solidago altissima",
                "Solidago canadensis",
                "Solidago speciosa",
                "Taraxacum officinale",
                "Tradescantia occidentalis",
                "Tragopogon dubius (major)",
                "Trifolium hybridum",
                "Trifolium pratense",
                "Trifolium repens",
                "Ulmus americana" , 
                "Verbascum thapsus", 
                "Viola pedata",
                "Viola sagittata",
                "Aster ericoides",
                "Bromus kalmii",
                "Chenopodium album",
                "Galium boreale",
                "Helianthemum bicknellii",
                "Helianthus giganteus",
                "Helianthus laetiflorus",
                "Hieracium aurantiacum",
                "Hieracium pilosum"  ,
                "Hypericum majus",
                "Lathyrus venosus",
                "Lechea stricta",
                "Lithospermum canescens",
                "Lithospermum caroliniense",
                "Panicum capillare",
                "Panicum praecocious",
                "Penstemon gracilis",
                "Polygala polygama",
                "Prunus virginiana",
                "Ranunculus rhomboideus",
                "Rhus glabra",
                "Rhus radicans",
                "Rudbeckia hirta",
                "Rudbeckia serotina (hirta)",
                "Setaria lutescens glauca",
                "Sisyrinchium angustifolium",
                "Sisyrinchium campestre",
                "Solidago missouriensis",
                "Toxicodendron radicans",
                "Viburnum lentago",                                   
                "Vicia villosa",
                "Viola pedatifida",
                "Zizia aurea")

sp_family <- c("Agrostis sp.",
               "Arabis sp.",
               "Aster sp.",
               "Bromus sp.",
               "Carex sp.",
               "Cyperus sp",
               "Cyperus sp.",
               "Digitaria sp.",
               "Elymus sp.",
               "Equisetum sp.",
               "Erigeron sp",
               "Festuca sp.",
               "Helianthus sp.",
               "Liatris sp.",
               "Lupinus sp.",
               "Oenothera sp.",
               "Oxalis sp.",
               "Oxybaphus sp.",
               "Panicum sp.",
               "Penstemon sp.",
               "Polygonatum sp.",
               "Plantago sp.",
               "Rhus sp.",
               "Rumex sp.",
               "Silene sp.",
               "Solidago sp.",
               "Sporobolus sp.",
               "Tradescantia sp.",
               "Trifolium sp.",
               "Anemone",
               "Aristida sp.",
               "Eragrostis sp.",
               "Hieracium",
               "Lactuca sp.",
               "Lithospermum sp.",
               "Lychia",
               "Prunus sp.",
               "Rosa",
               "Sisyrinchium",
               "Viola sp.")

sp_unid <- c("1st year woody",
             "32 species weeds",
             "Anro",
             "C3 grasses",
             "C4 grasses",
             "Cyperus tectorum", # probably Crepis tectorum?
             "Forbes",
             "Legumes",
             "Miscellaneous grasses",
             "Miscellaneous herbs",
             "Miscellaneous sp.",
             "Miscellaneous forbs",
             "Miscellaneous legumes",
             "Miscellaneous sp. 2",
             "Miscellaneous woody plants",
             "Sedges",
             "Swch",
             "Unknown sp.",
             "Unknown forb",
             "Unknown woody",
             "Unsorted biomass",
             "Miscellaneous sedge",
             "Real weeds",
             "Unknown forb 1",
             "Woody")

sp_noplant <- c("Bare ground",
                "Lichens",
                "Fungi",
                "Miscellaneous litter",
                "Miscellaneous woody litter",
                "Moss",
                "Mosses",
                "Mosses & lichens",
                "Mosses 2",
                "Mosses 3",
                "Rhus glabra litter",
                "Gopher mounds",
                "Woody debris",
                "Woody litter")



### Summed relative Biomass of Species with trait value ####
data_quality_Bio <- biom_ER_Bio %>%
  filter(Planted %in% "planted") %>%
  filter(!Mass.g.m.2. %in% c(0, NA, NaN)) %>%
  filter(!Year %in% c(2008, 2019, 2020)) %>%
  select(!c(X, Date, Planted)) %>%
  mutate(Species_Diaz = Species,
         Species_Diaz = gsub(Species_Diaz, pattern = "Achillea millefolium (lanulosa)", replacement = "Achillea millefolium"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Ambrosia coronopifolia",          replacement = "Ambrosia psilostachya"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rubus alleghaniensis",            replacement = "Rubus allegheniensis"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Setaria lutescens (glauca)",      replacement = "Cenchrus americanus"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Taraxacum officinale",            replacement = "Taraxacum campylodes"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Agropyron repens",                replacement = "Elymus repens"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Agropyron smithii",               replacement = "Elymus smithii"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Ambrosia artemisiifolia elatior", replacement = "Ambrosia artemisiifolia"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Aster azureus",                   replacement = "Symphyotrichum oolentangiense"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Baptisia leucantha",              replacement = "Baptisia alba"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Leptoloma cognatum",              replacement = "Digitaria cognata"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Lychnis alba",                    replacement = "Silene latifolia"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Panicum oligosanthes",            replacement = "Dichanthelium oligosanthes"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Petalostemum purpureum",          replacement = "Dalea purpurea"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Petalostemum villosum",           replacement = "Dalea villosa"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Polygonum convolvulus",           replacement = "Fallopia convolvulus"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rudbeckia serotina",              replacement = "Echinacea serotina"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Smilacina stellata",              replacement = "Maianthemum stellatum")) %>%
  merge(., 
        DiazTraits %>%
          select(species, Leaf_Area, Nmass, LMA, Plant_height, Diaspore_mass, LDMC),
        by.x = "Species_Diaz", by.y = "species", all.x = T) %>%
  # remove tree seedlings
  filter(!Species_Diaz %in% c("Quercus ellipsoidalis", "Quercus macrocarpa")) %>%
  group_by(BigBioPlot, PlantSpNum, Plot, Treatment, Year) %>%
  mutate(rel_biomass = Mass.g.m.2./sum(Mass.g.m.2.)) %>%
  summarize(
    # CWM_LeafArea     = weighted.mean(Leaf_Area,     Mass.g.m.2., na.rm = T),
    # CWM_Nmass        = weighted.mean(Nmass,         Mass.g.m.2., na.rm = T),
    # CWM_LMA          = weighted.mean(LMA,           Mass.g.m.2., na.rm = T),
    # CWM_PlantHeight  = weighted.mean(Plant_height,  Mass.g.m.2., na.rm = T),
    # CWM_DiasporeMass = weighted.mean(Diaspore_mass, Mass.g.m.2., na.rm = T),
    # CWM_LDMC         = weighted.mean(LDMC,          Mass.g.m.2., na.rm = T),
    
    rel_biomass_Nmass        = sum(rel_biomass[!is.na(Nmass)]),
    rel_biomass_LeafArea     = sum(rel_biomass[!is.na(Leaf_Area)]),
    rel_biomass_LMA          = sum(rel_biomass[!is.na(LMA)]),
    rel_biomass_LDMC         = sum(rel_biomass[!is.na(LDMC)]),
    rel_biomass_PlantHeight  = sum(rel_biomass[!is.na(Plant_height)]),
    rel_biomass_DiasporeMass = sum(rel_biomass[!is.na(Diaspore_mass)])
  ) %>%
  filter(!PlantSpNum %in% 1) %>%
  ungroup()

data_quality_OF <-  biom_ER_OF %>%
  filter(!Subplot %in% "West") %>%
  filter(!Mass.g.m.2. %in% c(0, NA, NaN)) %>%
  filter(!Year %in% c(2008, 2019, 2020)) %>%
  filter(!Species %in% c(sp_noplant)) %>%
  select(!c(X)) %>%
  mutate(Species_Diaz = Species,
         Species_Diaz = gsub(Species_Diaz, pattern = "Achillea millefolium (lanulosa)", replacement = "Achillea millefolium"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Ambrosia coronopifolia",          replacement = "Ambrosia psilostachya"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rubus alleghaniensis",            replacement = "Rubus allegheniensis"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Setaria lutescens (glauca)",      replacement = "Cenchrus americanus"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Taraxacum officinale",            replacement = "Taraxacum campylodes"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Agropyron repens",                replacement = "Elymus repens"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Agropyron smithii",               replacement = "Elymus smithii"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Ambrosia artemisiifolia elatior", replacement = "Ambrosia artemisiifolia"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Aster azureus",                   replacement = "Symphyotrichum oolentangiense"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Baptisia leucantha",              replacement = "Baptisia alba"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Leptoloma cognatum",              replacement = "Digitaria cognata"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Lychnis alba",                    replacement = "Silene latifolia"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Panicum oligosanthes",            replacement = "Dichanthelium oligosanthes"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Petalostemum purpureum",          replacement = "Dalea purpurea"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Petalostemum villosum",           replacement = "Dalea villosa"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Polygonum convolvulus",           replacement = "Fallopia convolvulus"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rudbeckia serotina",              replacement = "Echinacea serotina"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Smilacina stellata",              replacement = "Maianthemum stellatum"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Aster ericoides",                 replacement = "Symphyotrichum ericoides"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rhus radicans",                   replacement = "Toxicodendron radicans"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Rudbeckia serotina (hirta)",      replacement = "Echinacea serotina"),
         Species_Diaz = gsub(Species_Diaz, pattern = "Setaria lutescens glauca",        replacement = "Cenchrus americanus")) %>%
  merge(., 
        DiazTraits %>%
          select(species, Leaf_Area, Nmass, LMA, Plant_height, Diaspore_mass, LDMC),
        by.x = "Species_Diaz", by.y = "species", all.x = T) %>%
  # remove tree seedlings
  filter(!Species_Diaz %in% c("Quercus ellipsoidalis", "Quercus macrocarpa", "Prunus virginiana")) %>%
  group_by(Plot, Treatment, Year) %>%
  mutate(rel_biomass = Mass.g.m.2./sum(Mass.g.m.2.)) %>%
  summarize(
    # CWM_LeafArea     = weighted.mean(Leaf_Area,     Mass.g.m.2., na.rm = T),
    # CWM_Nmass        = weighted.mean(Nmass,         Mass.g.m.2., na.rm = T),
    # CWM_LMA          = weighted.mean(LMA,           Mass.g.m.2., na.rm = T),
    # CWM_PlantHeight  = weighted.mean(Plant_height,  Mass.g.m.2., na.rm = T),
    # CWM_DiasporeMass = weighted.mean(Diaspore_mass, Mass.g.m.2., na.rm = T),
    # CWM_LDMC         = weighted.mean(LDMC,          Mass.g.m.2., na.rm = T),
    
    rel_biomass_Nmass        = sum(rel_biomass[!is.na(Nmass)]),
    rel_biomass_LeafArea     = sum(rel_biomass[!is.na(Leaf_Area)]),
    rel_biomass_LMA          = sum(rel_biomass[!is.na(LMA)]),
    rel_biomass_LDMC         = sum(rel_biomass[!is.na(LDMC)]),
    rel_biomass_PlantHeight  = sum(rel_biomass[!is.na(Plant_height)]),
    rel_biomass_DiasporeMass = sum(rel_biomass[!is.na(Diaspore_mass)])
  ) %>%
  
  ungroup()


#___________________________________________________________________________####
## ASSEMBLE FIGURES FROM SUBFIGURES ####
### -> Figure 1 PCA ####
plot_grid(
  compplot + 
    theme(legend.position = c(0,0), 
          legend.justification = c(0,0), 
          legend.margin = margin(l = 10, b = 10)), # +
  # geom_rect(aes(xmin = ggplot_build(compplot_trajectory_details3)$layout$panel_params[[1]]$x.range[1],
  #               xmax = ggplot_build(compplot_trajectory_details3)$layout$panel_params[[1]]$x.range[2], 
  #               ymin = ggplot_build(compplot_trajectory_details3)$layout$panel_params[[1]]$y.range[1],
  #               ymax = ggplot_build(compplot_trajectory_details3)$layout$panel_params[[1]]$y.range[2],), 
  #           fill = NA, color = "black"),  
  compplot_trajectory_details3  +
    guides(alpha = guide_legend(order = 1),
           color = guide_legend(order = 2, nrow = 3)) +
    theme(legend.position = c(0,0),
          legend.justification = c(0,0),
          legend.margin = margin(l = 10, b = 10),
          legend.box = "center",
          legend.box.just = "left"),
  rel_heights = c(1,2), nrow = 2,
  axis = "rl",   align = "hv",
  labels = c("a)", "b)")
)

### -> Figure 2 Treatment effects 2024 ####
legend2 <- ggpubr::get_legend(
  plots_TrtEff_OF[["DiasporeMass"]]+ theme(legend.position = "bottom", legend.justification="center") + guides(linetype = "none")
)


plot_grid(
  plot_grid(
    
    plots_effect_correlations[["Nmass"]] + theme(legend.position = "none"),
    plots_effect_correlations[["DiasporeMass"]] + theme(legend.position = "none"),
    plots_effect_correlations[["PlantHeight"]] + theme(legend.position = "none"),
    plots_effect_correlations[["LeafArea"]] + theme(legend.position = "none"),
    plots_effect_correlations[["LMA"]] + theme(legend.position = "none"),
    plots_effect_correlations[["LDMC"]] + theme(legend.position = "none"),
    
    ncol=2, align= "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  ggdraw(legend2),
  ncol = 1,
  rel_heights = c(0.95, 0.05)
)

### -> Figure 3 Species biomass change####
heatmap + theme(legend.key.width = unit(1.2, "cm"))

### -> Figure S1 Detail PCA ~ Time####
compplot_trajectory_details

### -> Figure S2 Biodiversity CWM tissue N ####
treatment_legend <- ggpubr::get_legend(
  ggplot(CWM_OF, aes(y = CWM_LeafArea_Biom , x = Year, color = Treatment, fill = Treatment)) +
    theme(legend.position = "bottom", legend.justification="center") +
    scale_color_manual(breaks = mybreaks,
                       values = myvalues,
                       labels = mylabels) +
    scale_fill_manual(breaks = mybreaks,
                      values = myvalues,
                      labels = mylabels) +
    geom_line() 
)

year_legend <- ggpubr::get_legend(
  plots_Bio$DiasporeMass$SR_Year + 
    theme(legend.position = "bottom", legend.justification="center") +
    labs(color = "")
)

plot_grid(
  plot_grid(
    plots_Bio[["Nmass"]]$Treatment_Year + theme(legend.position = "none"),
    plots_FG_Bio[["Nmass"]]$Treatment_Year + theme(legend.position = "none"),
    plots_Bio[["Nmass"]]$Treatment_SR + theme(legend.position = "none"),   
    plots_FG_Bio[["Nmass"]]$Treatment_SR + theme(legend.position = "none"),
    plots_Bio[["Nmass"]]$SR_Year + theme(legend.position = "none"),        
    plots_FG_Bio[["Nmass"]]$SR_Year + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  plot_grid(ggdraw(treatment_legend), ggdraw(year_legend), rel_widths = c(3,1)),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S3 Biodiversity CWM leaf area ####
plot_grid(
  plot_grid(
    plots_Bio[["LeafArea"]]$Treatment_Year + theme(legend.position = "none"),
    plots_FG_Bio[["LeafArea"]]$Treatment_Year + theme(legend.position = "none"),
    plots_Bio[["LeafArea"]]$Treatment_SR + theme(legend.position = "none"),   
    plots_FG_Bio[["LeafArea"]]$Treatment_SR + theme(legend.position = "none"),
    plots_Bio[["LeafArea"]]$SR_Year + theme(legend.position = "none"),        
    plots_FG_Bio[["LeafArea"]]$SR_Year + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  plot_grid(ggdraw(treatment_legend), ggdraw(year_legend), rel_widths = c(3,1)),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S4 Biodiversity CWM plant height ####
plot_grid(
  plot_grid(
    plots_Bio[["PlantHeight"]]$Treatment_Year + theme(legend.position = "none"),
    plots_FG_Bio[["PlantHeight"]]$Treatment_Year + theme(legend.position = "none"),
    plots_Bio[["PlantHeight"]]$Treatment_SR + theme(legend.position = "none"),   
    plots_FG_Bio[["PlantHeight"]]$Treatment_SR + theme(legend.position = "none"),
    plots_Bio[["PlantHeight"]]$SR_Year + theme(legend.position = "none"),        
    plots_FG_Bio[["PlantHeight"]]$SR_Year + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  plot_grid(ggdraw(treatment_legend), ggdraw(year_legend), rel_widths = c(3,1)),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S5 Oldfield CWM plant height ####
plot_grid(
  plot_grid(
    plots_OF[["PlantHeight"]] + theme(legend.position = "none"),
    plots_FG_OF[["PlantHeight"]] + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)")),
  
  ggdraw(treatment_legend),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S6 Oldfield CWM LMA ####
plot_grid(
  plot_grid(
    plots_OF[["LMA"]] + theme(legend.position = "none"),
    plots_FG_OF[["LMA"]] + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)")),
  
  ggdraw(treatment_legend),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S7 Oldfield CWM LDMC ####
plot_grid(
  plot_grid(
    plots_OF[["LDMC"]] + theme(legend.position = "none"),
    plots_FG_OF[["LDMC"]] + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)")),
  
  ggdraw(treatment_legend),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S8 Treatment effects over time ####
plot_grid(
  
  plot_grid(
    plots_TrtEff[["Nmass"]] + theme(legend.position = "none"),
    plots_TrtEff[["DiasporeMass"]] + theme(legend.position = "none"),
    plots_TrtEff[["PlantHeight"]] + theme(legend.position = "none"),
    plots_TrtEff[["LeafArea"]] + theme(legend.position = "none"),
    plots_TrtEff[["LMA"]] + theme(legend.position = "none"),
    plots_TrtEff[["LDMC"]] + theme(legend.position = "none"),
    
    ncol=2, align= "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  ggdraw(legend2),
  ncol = 1,
  rel_heights = c(0.95, 0.05)
)

### -> Figure S9 Biodiversity inverse simpson ####
plot_grid(
  plot_grid(
    plots_Bio[["invsimpson"]]$Treatment_Year + theme(legend.position = "none"),
    plots_FG_Bio[["invsimpson"]]$Treatment_Year + theme(legend.position = "none"),
    plots_Bio[["invsimpson"]]$Treatment_SR + theme(legend.position = "none"),   
    plots_FG_Bio[["invsimpson"]]$Treatment_SR + theme(legend.position = "none"),
    plots_Bio[["invsimpson"]]$SR_Year + theme(legend.position = "none"),        
    plots_FG_Bio[["invsimpson"]]$SR_Year + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)", "c)", "d)", "e)", "f)")),
  
  plot_grid(ggdraw(treatment_legend), ggdraw(year_legend), rel_widths = c(3,1)),
  ncol = 1,
  rel_heights = c(0.9, 0.1))

### -> Figure S10 Oldfield inverse simpson ####
plot_grid(
  plot_grid(
    plots_OF[["invsimpson"]] + theme(legend.position = "none"),
    plots_FG_OF[["invsimpson"]] + theme(legend.position = "none"),
    
    ncol = 2, align = "hv",
    
    labels = c("a)", "b)")),
  
  ggdraw(treatment_legend),
  ncol = 1,
  rel_heights = c(0.9, 0.1))


### -> Figure S11 Species biomass response ~ trait correlation ####
plot_sp_TrtEff2

### -> Figure S12 Species biomass response correlation between experiments ####
plot_sp_TrtEff_cor


#___________________________________________________________________________####
## SAVE RDATA ####

save(list = ls(), 
     file = paste(EnemyRemovalDIR, "/results/CWMs_BIOM.RData", sep = ""))

