# 2025ConsumerLegacies
Code accompanying the manuscript “Legacy effects shape plant community responses to consumer exclusion”

The code data.assembly.R is used to assemble the data products from raw data that are analyzed in the analysis.R code. The text for the manuscipt is created in the MS_EnemyRemoval_CommunityShift.Rmd rmarkdown document. It sources the results that are produced by analysis.R.



## The following raw data and auxiliary files were used


### BigBio_All_Plots_info.csv
This file contains the design of the big biodiversity experiment at Cedar Creek
- Exp				Cedar Creek internal experiment number (120)
- Plot				Unique plot number within the big biodiversity experiment
- NumSp				Number of species planted in plot
- FgNum				Number of functional groups planted in plot
- Fgset				Code for unique species composition (?)
- C3				Whether (1) or not (0) C3 grasses were planted in plot
- C4				Whether (1) or not (0) C4 grasses were planted in plot
- Forb				Whether (1) or not (0) non-leguminous forbs were planted in plot
- Legume			Whether (1) or not (0) legumes were planted in plot
- Woody				Whether (1) or not (0) woody species were planted in plot
- Achillea millefolium(lanulosa)Whether (1) or not (0) the species Achillea millefolium was planted in plot
- ...				Whether (1) or not (0) the species ... was planted in plot
- Zizia aurea			Whether (1) or not (0) the species Zizia aurea was planted in plot


### E244Data_PlanFilev1.csv
This file contains the design of the consumer removal experiment nested within the big biodiversity experiment at Cedar Creek
- BigBioPlot		Big biodiversity experiment plot within which the experimental plot is nested. Corresponds to Plot in BigBio_All_Plots_Info.csv
- ER244.Plot		Unique plot number within the consumer removal experiment 
- PlantSpNum		Number of species planted in plot. Corresponds to NumSp in BigBio_All_Plots_Info.csv
- FieldColorCode	Plot marking color (per treatment)
- Insecticide		Whether (1) or not (0) insecticide is applied in the plot
- SoilDrenchFungicide	Whether (1) or not (0) soil drench fungicide is applied in the plot
- FoliarFungicide	Whether (1) or not (0) foliar fungicide is applied in the plot
- TreatmentName		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides)


### Cedar Creek Plant Taxon List.csv
This file contains information to all species entries across all Cedar Creek experiments
- Species			Species name entered in data sheets (incl. typos)
- Functional group		Abbreviated functional group abbreviation 
- Duration			Whether the plant is perennial or annual
- Lifeform			Life/growth form
- Pathway			Photosynthetic pathway
- Taxon				Unique taxon number
- Specid			Species name long abbreviation
- 5Lspecid			Species name short abbreviation
- Origin			Whether the plant is native or introduced
- ITISTaxon			ITIS Taxon number
- ITISRecognizedName		ITIS Recognized species name
- Family			Taxonomic family
- CountSpecies			-
- New name (USDA)		New species name
- Comment			Comments


### e244_Plant aboveground biomass data.txt
This file contains the (mostly) species level biomass data per plot from 2007-2019 from the enemy removal experiment in the big biodiversity experiment at Cedar Creek. It's header contains conditions for data re-use. 

Kinkel, Linda (2024): Plant aboveground biomass data: Natural Enemies, Plant Diversity and Plant Community Composition. DOI: 10.6073/pasta/c1ffd1c48b504e97a93fbe3583b2508f

- Year			Year when biomass was collected
- Date			Date when biomass was collected
- Plot			Unique plot number within the consumer removal experiment in the biodiversity experiment.
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides)
- Species		Species name
- Mass g/m^2)		Plant biomass in g/m^2


### e245_Plant aboveground biomass data.txt
This file contains the (mostly) species level biomass data per plot from 2007-2019 from the enemy removal experiment in the old agricultural field (oldfield) at Cedar Creek. It's header contains conditions for data re-use. 

Kinkel, Linda (2024): Plant aboveground biomass data: The influence of natural enemies on plant community composition and productivity. DOI: 10.6073/pasta/303607d5f92929a4b20ba127c47d21f0

- Year			Year when biomass was collected
- Date			Date when biomass was collected
- Plot			Unique plot number in the consumer removal experiment in the oldfield.	
- Subplot		In 2019 the plots were split into two subplot, one fertilized, one unfertilized. For this study, only unfertilized subplots were used
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- FertTrt		Whether (y) or not (N) the subplot was fertilized
- Species		Species name
- Mass g/m^2)		Plant biomass in g/m^2


### e244 e245 2020-2024 Abovegroujnd Biomass for SC.csv
This file contains the species level biomass data per plot from 2020-2024 from both enemy removal experiment
- Year			Year when biomass was collected
- Experiment		E244 (enemy removal in the biodiversity experiment), E245 (enemy removal in the oldfield experiment)
- Plot			Plot number within each experiment
- Subplot		In 2019 the plots were split into two subplot, one fertilized, one unfertilized. For this study, only unfertilized subplots were used (only for the oldfield experiment E245)
- BigBioPlot#		Big biodiversity experiment plot within which the experimental plot is nested. (only for biodiversity experiment E244)
- Block#		Block number (only for oldfield experiment E245)
- TreatmentCode		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- BigBioNumSp		Number of species planted in the plot (only for biodiversity experiment E244)
- FertTrt		Whether (y) or not (N) the subplot was fertilized (only for oldfield experiment E245)
- Species		Species name
- Mass g/m^2)		Plant biomass in g/m^2


### CDR_traits_SDaiz.csv
Cedar Creek species subset of plant trait data from Diaz et al. 2022 assembled by Jesús N. Pinto-Ledezma (https://github.com/jesusNPL/CDR_phyfunc)

Díaz, Sandra; Kattge, Jens; Cornelissen, Johannes H. C.; Wright, Ian J.; Lavorel, Sandra; Dray, Stéphane et al. (2022): The global spectrum of plant form and function: enhanced species-level trait dataset. In Scientific data 9 (1), p. 755. DOI: 10.1038/s41597-022-01774-9.

order	
- family		Taxonomic family 
- genus			Taxonomic genus
- species		Species name
- Leaf_Area		Leaf area [mm^2]
- Nmass			Leaf tissue nitrogen content [%]
- LMA			Leaf mass per area [g.m^-2]
- Plant_height		Plant height [m]
- Diaspore_mass		Seed mass [mg]
- SSD_observed		-
- LDMC			Leaf dry matter content [g/g]
- SSD_imputed		-
- SSD_combined		-




## The following aggregated data products are created by the data.assembly.R code


### speciesinfo_extended.csv
Cedar Creek Plant Taxon List.csv, but with missing taxa added.


### E244_biomass_clean.csv
cleaned version of e244_Plant aboveground biomass data.txt


### E245_biomass_clean.csv
cleaned version of e245_Plant aboveground biomass data.txt


### e244_CWMTraits.csv
CWM of traits calculated from E244_biomass_clean.csv and CDR_traits_SDaiz.csv
- Plot				Unique plot number within the consumer removal experiment in the biodiversity experiment.
- Treatment			Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides)
- Year				Year when biomass was collected
- CWM_LeafArea_Biom_sownSp	CWM of leaf area [mm^2]
- CWM_Nmass_Biom_sownSp		CWM of tissue nitrogen content [%]
- CWM_LMA_Biom_sownSp		CWM of leaf mass per area [g.m^-2]
- CWM_PlantHeight_Biom_sownSp	CWM of plant height [m]
- CWM_DiasporeMass_Biom_sownSp	CWM of seed mass [mg]
- CWM_LDMC_Biom_sownSp		CWM of leaf dry matter content [g/g]
- BigBioPlot			Big biodiversity experiment plot within which the experimental plot is nested. Corresponds to Plot in BigBio_All_Plots_Info.csv
- NumSp				Number of species planted in plot.
- FgNum				Number of functional groups planted in plot
- Fgset				Code for unique species composition (?)
- C3				Whether (1) or not (0) C3 grasses were planted in plot
- C4				Whether (1) or not (0) C4 grasses were planted in plot
- Forb				Whether (1) or not (0) non-leguminous forbs were planted in plot
- Legume			Whether (1) or not (0) legumes were planted in plot
- Woody				Whether (1) or not (0) woody species were planted in plot
				

### e245_CWMTraits.csv
CWM of traits calculated from E245_biomass_clean.csv and CDR_traits_SDaiz.csv
- Plot				Unique plot number in the consumer removal experiment in the oldfield.
- Treatment			Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- Subplot			In 2019 the plots were split into two subplot, one fertilized, one unfertilized. For this study, only unfertilized subplots were used
- Year				Year when biomass was collected
- CWM_LeafArea_Biom		CWM of leaf area [mm^2]
- CWM_Nmass_Biom		CWM of tissue nitrogen content [%]
- CWM_LMA_Biom			CWM of leaf mass per area [g.m^-2]
- CWM_PlantHeight_Biom		CWM of plant height [m]
- CWM_DiasporeMass_Biom		CWM of seed mass [mg]
- CWM_LDMC_Biom			WM of leaf dry matter content [g/g]


### e244_CWMFGAbundance.csv
Summed biomass of each functional group in E244_biomass_clean.csv 
- Year			Year when biomass was collected
- Plot			Unique plot number within the consumer removal experiment in the biodiversity experiment.
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides)
- Functional.group	Functional group
- Mass.g.m.2.		Biomass in g.m^-2


### e245_CWMFGAbundance.csv
Summed biomass of each functional group in E245_biomass_clean.csv 
- Plot			Unique plot number in the consumer removal experiment in the oldfield.
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- Subplot		In 2019 the plots were split into two subplot, one fertilized, one unfertilized. For this study, only unfertilized subplots were used
- Year			Year when biomass was collected
- Functional.group	Functional group
- Mass.g.m.2.		Biomass in g.m^-2


### E244_div.csv
Inverse Simpson index calculated from E244_biomass_clean.csv 
- Year			Year when biomass was collected
- Plot			Unique plot number within the consumer removal experiment in the biodiversity experiment.
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides)
- invsimpson_biomass	inverse Simpson index
- BigBioPlot		Big biodiversity experiment plot within which the experimental plot is nested. Corresponds to Plot in BigBio_All_Plots_Info.csv
- NumSp			Number of species planted in plot.


### E245_div.csv
Inverse Simpson index calculated from E245_biomass_clean.csv 
- Plot			Unique plot number in the consumer removal experiment in the oldfield.
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- Subplot		In 2019 the plots were split into two subplot, one fertilized, one unfertilized. For this study, only unfertilized subplots were used
- Year			Year when biomass was collected
- invsimpson_biomass	inverse Simpson index


### species_biomass_change.csv
biomass change of abundant species (occurred >10 times in any treatment across time) per treatment over time estimated with emtrends from lmer model Mass.g.m.2.~ Treatment * Year2 + (1|Plot) + (1|Year)
- Species		Species name
- Treatment		Unique treatment name (Control, Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- Year2			mean standardized Year
- Year2.trend		biomass change over time (slope) of species in treatment
- SE			standard error
- df			degrees of freedom
- lower.CL		lower end of 95% confidence interval
- upper.CL		upper end of 95% confidence interval 
- t.ratio		t ratio (is change over time different from 0)
- p.value		p value (is change over time different from 0)
- exp			experiment (Biodiversity, Oldfield)
- Functional.group	Functional group of species


### species_biomass_TrtEff.csv
treatment effect on biomass of abundant species (occurred >10 times in any treatment across time) per treatment over time estimated with emmeans from lmer model Mass.g.m.2.~ Treatment * Year2 + (1|Plot) + (1|Year)
- Species_Diaz		Species name in Diaz et al. 2022 data
- Species		Species name
- contrast		Compared treatments (Control vs. Treatment)
- estimate		difference between compared treatments
- SE			standard error
- df			degrees of freedom
- t.ratio		t ratio (does consumer removal change species biomass relative to the control)
- p.value		p value (does consumer removal change species biomass relative to the control)
- exp			experiment (Biodiversity, Oldfield)
- Functional.group	Functional group of species	
- Treatment		Unique treatment name (Insecticide, SoilDrenchFungicide, FoliarFungicide, AllPesticides, Fence)
- Leaf_Area		Leaf area [mm^2]
- Nmass			Leaf tissue nitrogen content [%]
- LMA			Leaf mass per area [g.m^-2]
- Plant_height		Plant height [m]
- Diaspore_mass		Seed mass [mg]
- LDMC			Leaf dry matter content [g/g]


