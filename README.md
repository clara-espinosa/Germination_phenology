README.md
================

[![DOI](https://zenodo.org/badge/440159075.svg)](https://zenodo.org/doi/10.5281/zenodo.11242675)

<figure>
<img src="microsites.png"
alt="Fellfield and snowbed microsites from the Cantabrian Mountains" />
<figcaption aria-hidden="true">Fellfield and snowbed microsites from the
Cantabrian Mountains</figcaption>
</figure>

# Microclimate regulates germination phenology in alpine plant communities

This repository stores all the information related to the manuscript
[*Microclimatic variation regulates seed germination phenology in alpine
plant communities*](add%20doi%20once%20manuscript%20is%20accepted),
including the raw datasets, the scripts to perform data cleaning,
analysis and figures.

## Contents

This repository is organised following the advice of [Wilson et
al. 2017](https://doi.org/10.1371/journal.pcbi.1005510) for recording
and storing research projects.

The following materials are available in the folders of this repository:

- `data` Data files including the [Raw germination data from the
  laboratory
  experiments](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/clean%20data.csv)
  and [Germination data from the field sowing
  experiments](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/field_germination.csv),
  [Species
  information](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/all_info.csv)
  datasets; [CHELSA bioclimatic variables for the study
  system](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/Bioclim_study_area.csv)
  for figure 1b, [Temperatures recorded in the field from 2008 to
  2019](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/temp_picos_raw.csv)
  then transformed into [Weekly temperature records for
  visualization](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/weekly_picos_graph.csv)
  in figure 1c and [Temperatures programs used in growth
  chambers](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/date_temp.csv)
  for figure 1d. It also includes the [Mean germination trait values
  table](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/meanvalues_graph.csv)
  and the [Efects sizes from MCMC-GLMM
  models](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/test_effectsize.csv)
  for figure 4. Finally a [Metadata
  file](https://github.com/clara-espinosa/Germination_phenology/blob/main/data/metadata.csv)
  explaining each csv variables.
- `doc` Files required for the submission at Ecology Letters,
- `results` Output of the `R` scripts, including figures,models and
  phylogenetic trees.
- `src` Scripts in `R` language used to clean the raw data, perform the
  analyses of the manuscript and create the figures.

## Abstract

- **Background** Understanding seed germination phenology is crucial for
  predicting plant responses to environmental changes. However, a
  substantial gap persists regarding how microclimatic conditions
  influence germination in seasonal ecosystems.  
- **Methods** We conducted a continuous seasonal experiment with fresh
  seeds to investigate germination phenology in 54 species from
  temperate and Mediterranean alpine communities. Using field
  microclimatic data series, we mimicked fellfield and snowbed
  conditions in growth chambers and we carried out field sowing
  experiments.  
- **Key Results** Both communities showed similar phenology responses to
  microclimatic variation, finding a consistent germination delay in
  snowbed compared to fellfield conditions. This effect was complemented
  by reduced dormancy and increased autumn germination in Mediterranean
  seeds.
- **Conclusions** Our results suggest a predictable phenological shift
  in the germination of alpine plants along microclimatic gradients. In
  warmer conditions with reduced snow cover, alpine species are expected
  to anticipate germination 52 days on average, with potential
  disrupting effects on cold-adapted plant communities.

## Citation

Please cite the repository, datasets and article as (to complete once
manuscript has been accepted for publication):
