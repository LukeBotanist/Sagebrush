# Sagebrush population genetics analysis pipeline

This repository hosts data and code for the population genetic analyses performed in Grossfurthner et al. (2023) Population structure and hybridization under contemporary and future climates in a heteroploid foundational shrub species (*Artemisia tridentata*). The free full text can be found [here](https://doi.org/10.3389/fpls.2023.1155868).

# Project description

This project aimed to identify the distinctiveness of *Artemisia tridentata* subspecies as well as the detection of hybridization. To this end, we sampled five transects throughout the western United States where a subspecies overlap was predicted using subspecies-specific climate niche models. Along each transect, we sampled multiple plots representing the parental and the potential hybrid habitats. We performed reduced representation sequencing and processed the data using a ploidy-informed genotyping approach. Here, I provide the basic codes used for the analyses of the heteroploid data.
Additional analyses, annotations and bug-fixes will follow soon. 

# Data analysis

All analyses are contained within R, however, the scripts provided includes links to external software (i.e. PolyRelatedness), which needs to be installed separately. A detailed description of each step can be found in the documentation directory `/docs/EBG_pipeline.pdf`.
The full pipeline `EBG_pipeline.R` can be executed using `Rscript EBG_pipeline.R` and will download data from this repository, setup directory structure and generate all tables and figures used for the analyses. 
**NOTE:** The software PolyRelatedness needs to be installed prior to analysis and the path to the PolyRelatedness directory needs to be set manually. Alternatively, to enable automatic download, extraction and execution of PolyRelatedness, corresponding lines in the code need to be uncommented.
This pipeline was written under R version 3.3.6 and tested up tp R version 4.3.1 on a aarch64-apple-darwin20 (64-bit) platform, running under: macOS Ventura 13.5.2.

- `data/` contains:
  - `EBG_data`: files with diploid and tetraploid genotypes, shared variants and Sample_ID
  - `PolyRelatedness_data`: required formatting files and example file
  - `STRUCTURE_data`: subspecies assignments and ancestry proportions as inferred in Grossfurthner et al. (2023)

- `docs/` contains tutorials for data analysis in different formats
  - `EBG_pipeline.html`
  - `EBG_pipeline.pdf`
  - `EBG_pipeline.md`
  - `EBG_pipeline_files`: directory with figures for renderind the markdown.

- `src/` contains:
  - `functions.R`: functions to convert and calculate data; invoked by the `EBG_pipeline.R script`
  - `EBG_pipeline.R`: script for data analysis and plotting; can be exectuted using `Rscript EBG_pipeline.R`. **Note:** The software PolyRelatedness needs to be installed first.



# Funding statement

Funding was provided by the NSF Idaho EPSCoR Program and by the National Science Foundation under award number OIA-1757324. Additional financial support was provided by the Stillinger Expedition fund, and the USDA Rocky Mountain Research Station.


# [Getting Started](docs/EBG_pipeline.md)
