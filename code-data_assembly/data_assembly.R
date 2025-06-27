# created 28.09.2023 by seraina.cappelli@gmail.com
# Project PI: Elizabeth Borer

# Kaiser Chiefs - Ruby
# Prince - Kiss
# Roaman - Because I'm Spiritual (The Kundalini Song)
# Natalie Merchant - Which Side Are You on?
# Loreley & Me - To the Sun

# project path:
EnemyRemovalDIR <- getwd() %>% paste(., "/", sep = "")

# R studio version , R Version R-4.3.1

library(tidyverse)
# dplyr     1.1.3     # readr     2.1.4
# forcats   1.0.0     # stringr   1.5.0
# ggplot2   3.4.3     # tibble    3.2.1
# lubridate 1.9.2     # tidyr     1.3.0
# purrr     1.0.2
library(vegan) # 2.6-4
library(lme4)
library(car)
library(emmeans)

#### GET DATA ####
##### experiment design and species information ####
plotinfo_BigBio <- read.csv("data-raw/BigBio_All_Plots_info.csv") %>%
  select(!c(Exp)) 
plotinfo_Bio <- read.csv("data-raw/E244Data_PlanFilev1.csv")
speciesinfo <- read.csv(paste(EnemyRemovalDIR, "/data-raw/Cedar Creek Plant Taxon List.csv", sep = "")) %>%
  mutate(Species = ifelse(Species %in% "Achillea millefolium(lanulosa)", "Achillea millefolium (lanulosa)", Species),
         Species = ifelse(Species %in% "Taraxicum officinalis", "Taraxacum officinale", Species),
         Species = ifelse(Species %in% "Rhus radicans", "Toxicodendron radicans", Species)) %>%
  mutate(Functional.group = as.factor(Functional.group),
         Duration         = as.factor(Duration),
         Lifeform         = as.factor(Lifeform),
         Pathway          = as.factor(Pathway),
         Origin           = as.factor(Origin)) %>%
  rbind.data.frame(., 
                   c(Species = "Oxybaphus sp.",
                     Functional.group = "F",
                     Duration = "Perennial",
                     Lifeform = "Non-legume forb",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = "Native",
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Nyctaginaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Plantago sp.",
                     Functional.group = "F",
                     Duration = "Perennial",
                     Lifeform = "Non-legume forb",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = "Introduced",
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Plantaginaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Elymus sp.",
                     Functional.group = "F",
                     Duration = NA,
                     Lifeform = "Grass",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Poaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "C3 grasses",
                     Functional.group = "C3",
                     Duration = NA,
                     Lifeform = "Grass",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Poaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = ""))  %>%
  rbind.data.frame(., 
                   c(Species = "C4 grasses",
                     Functional.group = "C3",
                     Duration = NA,
                     Lifeform = "Grass",
                     Pathway = "C4",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Poaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Hieracium pilosum",
                     Functional.group = "F",
                     Duration = NA,
                     Lifeform = "Non-legume forb",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Poaceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "1st year woody",
                     Functional.group = "W",
                     Duration = NA,
                     Lifeform = NA,
                     Pathway = NA,
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = NA,
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Forbes",
                     Functional.group = "F",
                     Duration = NA,
                     Lifeform = NA,
                     Pathway = NA,
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = NA,
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Legumes",
                     Functional.group = "L",
                     Duration = NA,
                     Lifeform = NA,
                     Pathway = NA,
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = NA,
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")) %>%
  rbind.data.frame(., 
                   c(Species = "Liatris sp",
                     Functional.group = "F",
                     Duration = "Perennial",
                     Lifeform = "Non-legume forb",
                     Pathway = "C3",
                     Taxon = NA,
                     Specid = NA,
                     X5LSpecid = NA,
                     Origin = NA,
                     ITISTaxon = NA,
                     ITISRecognizedName = NA,
                     Family = "Asteraceae",
                     CountSpecies = NA,
                     New.name..USDA. = NA,
                     Comment = "")                   ) 

write.csv(speciesinfo, paste(EnemyRemovalDIR, "/data-derived/speciesinfo_extended.csv", sep = ""))

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
         Species2 = gsub(Species2, pattern ="Achillea millefolium lanulosa ", replacement = "Achillea millefolium (lanulosa)")) %>%
  merge(plotinfo_Bio %>% 
          select("BigBioPlot", "ER244.Plot"),
        by.x = "Plot", by.y = "BigBioPlot") %>%
  rename(BigBioPlot = Plot,
         Plot = ER244.Plot) %>%
  select(BigBioPlot, Plot, Species, Species2, Species_Diaz) %>%
  mutate(Plot_Species_Diaz = paste(Plot, Species_Diaz),
         Plot_Species = paste(Plot, Species2))


##### species filters ####
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

sp_moss_lichen <- c("Moss",
                    "Mosses",
                    "Mosses & lichens",
                    "Mosses 2",
                    "Mosses 3",
                    "Lichens")

sp_moss <- c("Moss",
             "Mosses",
             "Mosses 2",
             "Mosses 3")

sp_lichen <- c("Lichens")

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



##### traits ####
DiazTraits <- read.csv("data-raw/CDR_traits_SDaiz.csv") 


###### biomass ####
biom_both_20_24 <- read.delim("data-raw/e244 e245 2020-2024 Abovegroujnd Biomass for SC.csv", sep = ",") %>%
  rename("Treatment" = "TreatmentCode")

biom_ER_Bio <- read.delim("data-raw/e244_Plant aboveground biomass data.txt", skip = 48) %>%
  select(!Treatment) %>%
  ## Some datasets have different treatment assigned than what they should have, 
  ## lets exclusively rely on the plotinformation data!
  merge(plotinfo_Bio %>% select(ER244.Plot, TreatmentName) %>% rename("Treatment" = "TreatmentName"),
        by.x = "Plot", by.y = "ER244.Plot") %>%
  rbind.data.frame(biom_both_20_24 %>% 
                     filter(Experiment %in% "E244") %>% 
                     select(!Treatment) %>%
                     merge(plotinfo_Bio %>% select(ER244.Plot, TreatmentName) %>% rename("Treatment" = "TreatmentName"),
                           by.x = "Plot", by.y = "ER244.Plot") %>%
                     mutate(Date = NA) %>%
                     select(Year, Date, Plot, Treatment, Species, Mass.g.m.2.)) %>%
  mutate(Year = as.integer(paste(Year)),
         Treatment = factor(Treatment, 
                            levels = c("Control", "Insecticide", 
                                       "FoliarFungicide", 
                                       "SoilDrenchFungicide", 
                                       "AllPesticides")),
         Species = str_squish(Species),
         Species = str_to_sentence(Species),
         Species = as.factor(Species)) %>%
  filter(!Year %in% NA) %>%
  filter(!Treatment %in% NA)


biom_ER_OF <- read.delim("data-raw/e245_Plant aboveground biomass data.txt", skip = 51) %>%
  rbind.data.frame(biom_both_20_24 %>% 
                     filter(Experiment %in% "E245") %>% 
                     mutate(Date = NA,
                            Strip = NA) %>%
                     select(Year, Date, Plot, Subplot, Treatment, FertTrt, Species, Mass.g.m.2., Strip)) %>%
  mutate(Year = as.integer(paste(Year)),
         Plot = as.integer(paste(Plot)),
         Treatment = factor(Treatment, 
                            levels = c("Control", "Insecticide", 
                                       "FoliarFungicide", 
                                       "SoilDrenchFungicide", 
                                       "AllPesticides", "Fenced")), 
         FertTrt = factor(FertTrt),
         Species = factor(Species),
         Species = str_squish(Species),
         Species = str_to_sentence(Species),
         Species = as.factor(Species)) %>%
  filter(!Year %in% c(NA,1,2,3,4,5,6)) %>%
  # in 2017 and 2018 biomass was collected in two strips and values need to be 
  # averaged to get one measurement per plant and year
  group_by(Year, Plot, Subplot, Treatment, FertTrt, Species) %>%
  reframe(Mass.g.m.2. = ifelse(Year %in% c(2017, 2018), sum(Mass.g.m.2.)/2, sum(Mass.g.m.2.)))




#### CLEAN BIOMASS DATA ####
##### Biodiversity ####
biom_ER_Bio <- biom_ER_Bio %>%
  mutate(Species = recode_factor(Species,
                                 "Cypeus sp"                      = "Cyperus sp." ,
                                 "Achillea millefolium(lanulosa)" = "Achillea millefolium (lanulosa)",
                                 "Achillea millefolium (lanulosa) " = "Achillea millefolium (lanulosa)",
                                 "Miscellaneous forb"             = "Miscellaneous forbs",
                                 "Miscellaneous grass"            = "Miscellaneous grasses",
                                 "Taraxicum officinalis"          = "Taraxacum officinale")) %>%
  # some entries are double (identical biomass numbers). Most likely the same
  # biomass bags were weighted twice -> remove duplicates
  distinct(.keep_all = T) %>%
  # some entries are doubled, but without identical biomass number. Most likely
  # the species was accidentally sorted into two piles and not combined leater
  group_by(Year, Date, Plot, Treatment, Species) %>%
  summarize(Mass.g.m.2. = sum(Mass.g.m.2.))


##### Oldfield ####
biom_ER_OF <- biom_ER_OF %>%
  mutate(Species = recode_factor(Species,
                                 "Achillea millefolium(lanulosa)" = "Achillea millefolium (lanulosa)",
                                 "Achillea millefolium (lanulosa) " = "Achillea millefolium (lanulosa)",
                                 "Miscellaneous forb"             = "Miscellaneous forbs",
                                 "Miscellaneous grass"            = "Miscellaneous grasses",
                                 "Taraxicum officinalis"          = "Taraxacum officinale",
                                 "Viola petatifida"               = "Viola pedatifida",
                                 "Horsetail"                      = "Equisetum sp.", 
                                 "Rudbeckia serotina (hirta)"     = "Rudbeckia serotina",
                                 "Setaria glauca"                 = "Setaria lutescens (glauca)")) %>%
  # some entries are double (identical biomass numbers). Most likely the same 
  # biomass bags were weighted twice -> remove duplicates
  distinct(.keep_all = T)




write.csv(biom_ER_Bio %>%
            mutate(Plot_Species = paste(Plot, Species),
                   Planted = case_when(Plot_Species %in% design$Plot_Species ~ "planted", .default = "weed")) %>%
            merge(plotinfo_Bio %>% select(BigBioPlot, ER244.Plot, PlantSpNum), by.x = "Plot", by.y = "ER244.Plot") %>%
            select(Year, Date, BigBioPlot, PlantSpNum, Plot, Treatment, Species, Mass.g.m.2., Planted)
          , "data-derived/E244_biomass_clean.csv")
write.csv(biom_ER_OF, "data-derived/E245_biomass_clean.csv")


#### CALCULATE CWM TRAITS ####
##### Biodiversity ####
CWMTraits_ER_Bio_expSp_biomass <- biom_ER_Bio %>%
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
  mutate(Plot_Species_Diaz = paste(Plot, Species_Diaz)) %>%
  filter(Plot_Species_Diaz %in% design$Plot_Species_Diaz) %>%
  merge(., 
        DiazTraits %>%
          select(species, Leaf_Area, Nmass, LMA, Plant_height, Diaspore_mass, LDMC),
        by.x = "Species_Diaz", by.y = "species") %>%
  # remove tree seedlings
  filter(!Species_Diaz %in% c("Quercus ellipsoidalis", "Quercus macrocarpa")) %>%
  group_by(Plot, Treatment, Year) %>%
  summarize(CWM_LeafArea     = weighted.mean(Leaf_Area,     Mass.g.m.2., na.rm = T),
            CWM_Nmass        = weighted.mean(Nmass,         Mass.g.m.2., na.rm = T),
            CWM_LMA          = weighted.mean(LMA,           Mass.g.m.2., na.rm = T),
            CWM_PlantHeight  = weighted.mean(Plant_height,  Mass.g.m.2., na.rm = T),
            CWM_DiasporeMass = weighted.mean(Diaspore_mass, Mass.g.m.2., na.rm = T),
            CWM_LDMC         = weighted.mean(LDMC,          Mass.g.m.2., na.rm = T)) %>%
  merge(.,
        plotinfo_Bio %>%
          merge(plotinfo_BigBio %>% select(Plot, NumSp, FgNum, Fgset, C3, C4, Forb, Legume, Woody), by.x = "BigBioPlot", by.y = "Plot") %>%
          select(ER244.Plot, BigBioPlot, NumSp, FgNum, Fgset, C3, C4, Forb, Legume, Woody),
        by.x = "Plot",
        by.y = "ER244.Plot")



##### Oldfield #####
CWMTraits_ER_OF_Biom <- biom_ER_OF %>%
  
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
        by.x = "Species_Diaz", by.y = "species") %>%
  # remove trees
  filter(!Species_Diaz %in% "Prunus virginiana") %>%
  group_by(Plot, Treatment, Subplot, Year) %>%
  summarize(CWM_LeafArea     = weighted.mean(Leaf_Area,     Mass.g.m.2., na.rm = T),
            CWM_Nmass        = weighted.mean(Nmass,         Mass.g.m.2., na.rm = T),
            CWM_LMA          = weighted.mean(LMA,           Mass.g.m.2., na.rm = T),
            CWM_PlantHeight  = weighted.mean(Plant_height,  Mass.g.m.2., na.rm = T),
            CWM_DiasporeMass = weighted.mean(Diaspore_mass, Mass.g.m.2., na.rm = T),
            CWM_LDMC         = weighted.mean(LDMC,          Mass.g.m.2., na.rm = T))




write.csv(CWMTraits_ER_Bio_expSp_biomass%>%
            # note clumsy rename, because I tried a few different ways of calculating
            # CWM. This names are required for future code to run, even though, the other
            # CWM calculations are removed from the data assembly at this stage.
            rename(CWM_LeafArea_Biom_sownSp = CWM_LeafArea,
                   CWM_Nmass_Biom_sownSp = CWM_Nmass,
                   CWM_LMA_Biom_sownSp = CWM_LMA,
                   CWM_PlantHeight_Biom_sownSp = CWM_PlantHeight ,
                   CWM_LDMC_Biom_sownSp = CWM_LDMC,
                   CWM_DiasporeMass_Biom_sownSp  = CWM_DiasporeMass), 
          "data-derived/e244_CWMTraits.csv")

write.csv(CWMTraits_ER_OF_Biom %>%
            mutate(Subplot = case_when(Subplot %in% "na" ~ "", .default=Subplot),
                   Subplot = factor(Subplot)) %>%
            # note clumsy rename, because I tried a few different ways of calculating
            # CWM. This names are required for future code to run, even though, the other
            # CWM calculations are removed from the data assembly at this stage.
            rename(CWM_LeafArea_Biom = CWM_LeafArea,
                   CWM_Nmass_Biom = CWM_Nmass,
                   CWM_LMA_Biom = CWM_LMA,
                   CWM_PlantHeight_Biom = CWM_PlantHeight ,
                   CWM_LDMC_Biom = CWM_LDMC,
                   CWM_DiasporeMass_Biom  = CWM_DiasporeMass),
          "data-derived/e245_CWMTraits.csv")



# wrangle the plotinfo of Big Bio into a community matrix
sown_community_matrix_BigBio <- plotinfo_BigBio %>% 
  select(!c(NumSp, FgNum, Fgset, C3, C4, Forb, Legume, Woody)) %>%
  # remove columns with species for which we do not have R* values
  select(!c('Agropyron.smithii', 
            'Astragalus.canadensis', 
            'Baptisia.leucantha', 
            'Bare.ground', 
            'Bouteloua.gracilis', 
            'Bromus.inermis', 
            'Buchloe.dactyloides', 
            'Calamagrostis.canadensis', 
            'Elymus.canadensis', 
            'Leersia.oryzoides', 
            'Monarda.fistulosa',             
            'Quercus.ellipsoidalis', 
            'Quercus.macrocarpa', 
            'Sporobolus.cryptandrus', 
            'Trifolium.hybridum', 
            'Vicia.villosa',
            'Zizia.aurea')) %>%
  rename("Achillea millefolium (lanulosa)" = "Achillea.millefolium.lanulosa.") %>%
  rename_with(~ gsub(".", " ", .x, fixed = T)) %>%
  # add zero columns for species in the Rstar distance matrix that were not sown
  # in e120. Distance matrices are scaled based on max and min distance between
  # species. To make sown MPD R* comparable to measured MPD R* in both experi-
  # ments, I use the full distance matrix, even though, some of the species are
  # not part of the experiment. For the MPD function to run smooth, the distance
  # matrix needs to contain the exact same species as the community matrix.
  mutate('Agastache foeniculum'            = 0,
         'Agropyron repens'                = 0,
         'Agrostis scabra'                 = 0,
         'Ambrosia artemisiifolia elatior' = 0,
         'Amorpha canescens'               = 0,
         'Anemone cylindrica'              = 0,
         'Asclepias syriaca'               = 0,
         'Asclepias verticillata'          = 0,
         'Aster azureus'                   = 0,
         'Aster ericoides'                 = 0,
         'Bouteloua curtipendula'          = 0,
         'Calamovilfa longifolia'          = 0,
         'Coreopsis palmata'               = 0,
         'Desmodium canadense'             = 0,
         'Penstemon grandiflorus'          = 0,
         'Petalostemum candidum'           = 0,
         'Petalostemum villosum'           = 0,
         'Potentilla arguta'               = 0,
         'Rudbeckia serotina'              = 0,
         'Solidago nemoralis'              = 0,
         'Solidago speciosa'               = 0,
         'Solidago rigida'                 = 0,
         'Stipa spartea'                   = 0,
         Plot = paste("Plot", Plot, sep = "_")) %>%
  column_to_rownames("Plot") %>%
  as.matrix()





# FUNCTIONAL GROUP ABUNDANCES ####
##### Biodiversity ####
functgroup_biom_ER_Bio <- biom_ER_Bio %>% 
  mutate(Plot_Species = paste(Plot, Species)) %>%
  filter(Plot_Species %in% design$Plot_Species) %>%
  merge(speciesinfo %>% select(Species, Functional.group)) %>%
  group_by(Year, Plot, Treatment, Functional.group) %>%
  summarize(Mass.g.m.2.  = sum(Mass.g.m.2., na.rm = T ))  %>%
  
  rbind.data.frame(
    biom_ER_Bio %>%
      filter(Species %in% sp_moss_lichen) %>%
      mutate(Functional.group = "moss_lichen") %>%
      group_by(Year, Plot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T )) 
  )%>%
  
  rbind.data.frame(
    biom_ER_Bio %>%
      filter(Species %in% sp_moss) %>%
      mutate(Functional.group = "moss") %>%
      group_by(Year, Plot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T)) 
  )%>%
  
  rbind.data.frame(
    biom_ER_Bio %>%
      filter(Species %in% sp_lichen) %>%
      mutate(Functional.group = "lichen") %>%
      group_by(Year, Plot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T)) 
  )  %>% 
  
  # add zeros
  pivot_wider(
    values_from = Mass.g.m.2.,
    names_from = Functional.group,
    values_fill = 0
  ) %>% 
  
  pivot_longer(cols = c(C3, 'F', L, C4, W, 
                        moss_lichen,  moss),
               names_to = "Functional.group", 
               values_to = "Mass.g.m.2.")


##### Oldfield ####
functgroup_biom_ER_OF <- biom_ER_OF %>%
  merge(speciesinfo %>% select(Species, Functional.group)) %>%
  group_by(Year, Plot, Subplot, Treatment, Functional.group) %>%
  summarize(Mass.g.m.2.  = sum(Mass.g.m.2. )) %>%
  
  rbind.data.frame(
    biom_ER_OF %>%
      filter(Species %in% sp_moss_lichen) %>%
      mutate(Functional.group = "moss_lichen") %>%
      group_by(Year, Plot, Subplot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T )) 
  )%>%
  
  rbind.data.frame(
    biom_ER_OF %>%
      filter(Species %in% sp_moss) %>%
      mutate(Functional.group = "moss") %>%
      group_by(Year, Plot, Subplot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T)) 
  )%>%
  
  rbind.data.frame(
    biom_ER_OF %>%
      filter(Species %in% sp_lichen) %>%
      mutate(Functional.group = "lichen") %>%
      group_by(Year, Plot, Subplot, Treatment, Functional.group) %>%
      summarize(Mass.g.m.2. = sum(Mass.g.m.2., na.rm = T)) 
  ) %>% 
  
  # add zeros
  pivot_wider(
    values_from = Mass.g.m.2.,
    names_from = Functional.group,
    values_fill = 0
  ) %>% 
  
  pivot_longer(cols = c(Unknown,  C3, 'F', L, S, `NA`,C4, O, W, G, 
                        moss_lichen,  moss,  lichen),
               names_to = "Functional.group", 
               values_to = "Mass.g.m.2.")



write.csv(functgroup_biom_ER_Bio, "data-derived/e244_CWMFGAbundance.csv")
write.csv(functgroup_biom_ER_OF, "data-derived/e245_CWMFGAbundance.csv")


#### INVERSE SIMPSON ####
##### Biodiversity ####
community_matrix_ER_Bio_biom <- biom_ER_Bio %>%
  ungroup() %>%
  filter(!Species %in% c(sp_noplant)) %>%
  mutate(Plot_Species = paste(Plot, Species)) %>%
  filter(Plot_Species %in% design$Plot_Species) %>%
  mutate(Plot_year = paste(Plot, Year, Treatment, sep = "_")) %>%
  select(Plot_year, Species, Mass.g.m.2. )  %>%
  pivot_wider(id_cols     = Plot_year,
              names_from  = Species,
              values_from = Mass.g.m.2. ,
              values_fill = 0) %>%
  column_to_rownames("Plot_year") %>%
  select(sort(colnames(.)))

div_ER_Bio_biom <- community_matrix_ER_Bio_biom %>%
  diversity(index = "invsimpson") %>%
  data.frame() %>%
  rename("invsimpson_biomass" = ".") %>%
  rownames_to_column("plot_year_treatment") %>%
  separate_wider_delim(cols = plot_year_treatment,
                       delim = "_",
                       names = c("Plot", "Year", "Treatment")) %>%
  mutate(Plot = as.integer(Plot),
         Year = as.integer(Year),
         Treatment = factor(Treatment, levels = c("Control", 
                                                  "Fenced", 
                                                  "Insecticide", 
                                                  "SoilDrenchFungicide", 
                                                  "FoliarFungicide",
                                                  "AllPesticides"))) %>%
  merge(plotinfo_Bio %>% rename("NumSp"="PlantSpNum") %>% select(BigBioPlot, ER244.Plot, NumSp), by.x = "Plot", by.y = "ER244.Plot")


##### Oldfield ####
community_matrix_ER_OF_biom <- biom_ER_OF %>%
  filter(!Species %in% c(sp_noplant, "Miscellaneous litter", "Miscellaneous grasses",
                         "Miscellaneous herbs", "Fungi", "Mosses & lichens",
                         "Unknown forb 1", "Mosses", "Moss", "")) %>%
  mutate(Plot_year = paste(Plot, Year, Subplot, Treatment, sep = "_")) %>%
  select(Plot_year, Species, Mass.g.m.2.)  %>%
  pivot_wider(id_cols     = Plot_year,
              names_from  = Species,
              values_from = Mass.g.m.2.,
              values_fill = 0) %>%
  column_to_rownames("Plot_year")


div_ER_OF_biom <- community_matrix_ER_OF_biom %>%
  diversity(index = "invsimpson") %>%
  data.frame %>%
  rename("invsimpson_biomass" = ".") %>%
  rownames_to_column("plot_year_treatment") %>%
  separate_wider_delim(cols = plot_year_treatment,
                       delim = "_",
                       names = c("Plot", "Year", "Subplot", "Treatment")) %>%
  mutate(Plot = as.integer(Plot),
         Year = as.integer(Year),
         Subplot = as.factor(Subplot),
         Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides")))

write.csv(div_ER_Bio_biom, "data-derived/E244_div.csv")
write.csv(div_ER_OF_biom, "data-derived/E245_div.csv")




#### SPECIES BIOMASS CHANGE ####

# find how often a species occurs minimum per treatment across all plots and time points
# find how often a species occurs minimum per treatment across all plots and time points
abundant_species_per_treatment_Bio <- biom_ER_Bio %>%
  merge(plotinfo_Bio %>% select(ER244.Plot, PlantSpNum),
        by.x = "Plot", by.y = "ER244.Plot") %>%
  filter(!PlantSpNum %in% 1) %>%
  select(!PlantSpNum) %>%
  filter(!Mass.g.m.2. %in% c(0,NA)) %>%
  group_by(Species, Treatment) %>%
  summarize(N = length(Mass.g.m.2.)) %>%
  ungroup() %>%
  filter(N > 9.5) %>%
  filter(!Species %in% c(sp_noplant, sp_unid)) %>%
  mutate(sp_trt = paste(Species, Treatment))


abundant_species_per_treatment_OF <- biom_ER_OF %>%
  filter(!Subplot %in% "West") %>%
  filter(!Mass.g.m.2. %in% c(0, NA)) %>%
  group_by(Species, Treatment) %>%
  summarize(N = length(Mass.g.m.2.)) %>%
  ungroup() %>%
  filter(N > 9.5) %>%
  filter(!Species %in% c(sp_noplant, sp_unid)) %>%
  mutate(sp_trt = paste(Species, Treatment))



# add zeros to the data where necessary
biomass_for_sp_model_bio <- biom_ER_Bio %>%
  filter(Species %in% sp_species) %>%
  filter(!Year %in% c(2007, 2008, 2019, 2020)) %>%
  filter(!Mass.g.m.2. %in% NA) %>%
  mutate(sp_trt = paste(Species, Treatment)) %>%
  filter(sp_trt %in% abundant_species_per_treatment_Bio$sp_trt) %>%
  # select only those species that are sown in a plot
  mutate(Plot_Species = paste(Plot, Species)) %>%
  filter(Plot_Species %in% design$Plot_Species) %>%
  pivot_wider(id_cols     = c( "Plot", "Treatment", "Species"),
              names_from  = Year,
              values_from = Mass.g.m.2.,
              values_fn   = mean,
              values_fill = 0) %>%
  ungroup() %>%
  data.frame() %>%
  mutate(N_check = rowSums(select(.,-Plot, -Treatment, -Species))) %>%
  filter(N_check > 0) %>%
  pivot_longer(cols      = !c("Plot", "Treatment", "Species", "N_check"),
               names_to  = "Year",
               values_to = "Mass.g.m.2.") %>%
  mutate(Year = as.numeric(paste(gsub(Year, pattern = "X", replacement = "")))) %>%
  droplevels()



biomass_for_sp_model_OF <- biom_ER_OF %>%
  filter(Species %in% sp_species) %>%
  filter(!Subplot %in% "West") %>%
  filter(!Mass.g.m.2. %in% NA) %>%
  filter(!Year %in% c(2007, 2008, 2019, 2020)) %>%
  mutate(sp_trt = paste(Species, Treatment)) %>%
  filter(sp_trt %in% abundant_species_per_treatment_OF$sp_trt) %>%
  pivot_wider(id_cols     = c("Plot", "Treatment", "Species"),
              names_from  = Year,
              values_from = Mass.g.m.2.,
              values_fn   = mean,
              values_fill = 0) %>%
  ungroup() %>%
  data.frame() %>%
  mutate(N_check = rowSums(select(.,-Plot, -Treatment, -Species))) %>%
  filter(N_check > 0) %>%
  pivot_longer(cols      = !c("Plot", "Treatment", "Species", "N_check"),
               names_to  = "Year",
               values_to = "Mass.g.m.2.") %>%
  mutate(Year = as.numeric(paste(gsub(Year, pattern = "X", replacement = "")))) %>%
  droplevels()


# run models
posthoc_biom_slopes_bio <- data.frame()
posthoc_biom_TrtEff_bio <- data.frame() 
for (i in levels(droplevels(biomass_for_sp_model_bio$Species))){
  df = biomass_for_sp_model_bio %>%
    filter(Species %in% i) %>%
    mutate(Mass.g.m.2. = scale(Mass.g.m.2.))
  df$Year2 <- df$Year - 2007
  
  if (length(unique(df$Treatment)) > 1) {
    model = lmer(Mass.g.m.2.~ Treatment * Year2 + (1|Plot) + (1|Year), data = df)
    
    slopes = emtrends(model, ~ Treatment | Year2, var = "Year2", infer = T) %>%
      data.frame() %>%
      mutate(Species = paste(i))
    
    if ("Control" %in% df$Treatment) {
      TrtEff = (emmeans(model, pairwise~Treatment) )$contrasts %>%
        data.frame() %>%
        filter(grepl("Control", contrast)) %>%
        mutate(estimate = estimate * (-1),
               Species = paste(i))
    } else {
      TrtEff = data.frame(contrast = NA,
                          estimate = NA,
                          SE = NA,
                          df = NA,
                          t.ratio = NA,
                          p.value = NA,
                          Species = paste(i))
    }
    
  } else {
    model <- lmer(Mass.g.m.2. ~ Year2 + (1 | Plot) + (1 | Year), data = df)
    Treatment = unique(df$Treatment)
    
    slopes <- emtrends(model, ~ Year2, var = "Year2", infer = TRUE) %>%
      data.frame() %>%
      mutate(Species = i, Treatment = Treatment)
    
    TrtEff = data.frame(contrast = NA,
                        estimate = NA,
                        SE = NA,
                        df = NA,
                        t.ratio = NA,
                        p.value = NA,
                        Species = paste(i))
  }
  
  posthoc_biom_slopes_bio = rbind.data.frame(posthoc_biom_slopes_bio, slopes)
  posthoc_biom_TrtEff_bio = rbind.data.frame(posthoc_biom_TrtEff_bio, TrtEff)
}


posthoc_biom_slopes_of <- data.frame()
posthoc_biom_TrtEff_of <- data.frame()
for (i in levels(droplevels(biomass_for_sp_model_OF$Species))){
  df = biomass_for_sp_model_OF %>% 
    filter(Species %in% i) %>%
    mutate(Mass.g.m.2. = scale(Mass.g.m.2.))
  df$Year2 <- df$Year - 2007
  
  if (length(unique(df$Treatment)) > 1) {
    model = lmer(Mass.g.m.2. ~ Treatment * Year2 + (1|Plot) + (1|Year), data = df)
    
    slopes = emtrends(model, ~ Treatment | Year2, var = "Year2", infer = T) %>%
      data.frame() %>%
      mutate(Species = paste(i))
    
    if ("Control" %in% df$Treatment) {
      TrtEff = (emmeans(model, pairwise~Treatment) )$contrasts %>%
        data.frame() %>%
        filter(grepl("Control", contrast)) %>%
        mutate(estimate = estimate * (-1),
               Species = paste(i))
    } else {
      TrtEff = data.frame(contrast = NA,
                          estimate = NA,
                          SE = NA,
                          df = NA,
                          t.ratio = NA,
                          p.value = NA,
                          Species = paste(i))
    }
    
  } else {
    model <- lmer(Mass.g.m.2. ~ Year2 + (1 | Plot) + (1 | Year), data = df)
    
    Treatment = unique(df$Treatment)
    
    slopes <- emtrends(model, ~ Year2, var = "Year2", infer = TRUE) %>%
      data.frame() %>%
      mutate(Species = i, Treatment = Treatment)
    
    TrtEff = data.frame(contrast = NA,
                        estimate = NA,
                        SE = NA,
                        df = NA,
                        t.ratio = NA,
                        p.value = NA,
                        Species = paste(i))
  }
  
  posthoc_biom_slopes_of = rbind.data.frame(posthoc_biom_slopes_of, slopes)
  posthoc_biom_TrtEff_of = rbind.data.frame(posthoc_biom_TrtEff_of, TrtEff)
}


# merge results 

posthoc_slopes_biom <- 
  rbind.data.frame(
    posthoc_biom_slopes_bio %>%
      mutate(exp = "Biodiversity experiment"), 
    posthoc_biom_slopes_of %>% 
      mutate(exp = "Oldfield experiment")
  ) %>%
  
  merge(speciesinfo %>% 
          select(Species, Functional.group),
        by = "Species", 
        all.x = T) %>%
  arrange(Functional.group, Species) %>%
  mutate(Species = factor(Species, unique(Species)),
         Treatment = factor(Treatment, levels = c("Control", "Fenced", "Insecticide", "SoilDrenchFungicide", "FoliarFungicide", "AllPesticides")))


posthoc_TrtEff_biom <- 
  rbind.data.frame(
    posthoc_biom_TrtEff_bio %>%
      mutate(exp = "Biodiversity experiment"),
    posthoc_biom_TrtEff_of %>%
      mutate(exp = "Oldfield experiment")) %>%
  
  merge(speciesinfo %>% 
          select(Species, Functional.group),
        by = "Species", 
        all.x = T) %>%
  mutate(Species = factor(Species, unique(Species)),
         Treatment = gsub(contrast, pattern = "Control - ", replacement = ""),
         Treatment = factor(Treatment, levels = c("Control", 
                                                  "Fenced", 
                                                  "Insecticide", 
                                                  "SoilDrenchFungicide", 
                                                  "FoliarFungicide", 
                                                  "AllPesticides"))) %>%
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
        by.x = "Species_Diaz", by.y = "species", all.x =T)



write.csv(posthoc_slopes_biom, "data-derived/species_biomass_change.csv")
write.csv(posthoc_TrtEff_biom, "data-derived/species_biomass_TrtEff.csv")
