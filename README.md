# tufa_moss_diversity

## Repository Overview

This repository contains R scripts and processed data used for the manuscript:

**“Patterns and drivers of tufa moss diversity across alpine petrifying springs landscape gradient on the Tibetan Plateau”** (in review)

---

## Repository Structure

- **`env.csv`** contains the primary raw dataset used for all statistical analyses and visualizations. Additional details can be found in **`3_Supporting Table.xlsx`**, Sheet2 ("Variables data"), including:

  - Plot identifiers (`Sites`, `Plots`)
  - Geographic coordinates (`Lng`, `Lat`)
  - Habitat environmental variables (`Alt`, `Temp`, `RH`, `pH`, `Cond`, `OC`, `TN`, `TP`)
  - Moss morphological and physiological attributes (`MTC`, `MTN`, `MTP`, `Chl_a`, `Chl_b`, `Chl_a.b`, `Chl_a.b.1`, `MDA`, `TSS`, `POD`, `SOD`, `CAT`)
  - Diversity metrics (`Richness`, `Shannon`, `Simpson`, `NMDS1`)
  - Derived principal component variables (`PC1_diversity`, `PC1_geospatial`, `PC1_meteorological`, `PC1_Water`, `PC1_substrate`, `PC1_moss`, `PC1_photosynthesis`, `PC1_osmotic`, `PC1_enzyme`)

  The proportion of variance explained by the first principal component and the correlations between each principal component and its original variables are provided in **`3_Supporting Table.xlsx`**, Sheet4 ("PCA result").

  Detailed definitions of all variables, abbreviations, and descriptions are available in **`3_Supporting Table.xlsx`**, Sheet3 ("Data description").

- **`moss.csv`** contains moss community abundance data (cover values used as abundance proxies). Additional details are available in **`3_Supporting Table.xlsx`**, Sheet1 ("Tufa moss data"). This dataset is primarily used for community composition analyses and beta diversity assessments.

- **`3_Supporting Table.xlsx`**  
  - Sheet1–3 provide moss community abundance data, variable data, and detailed data descriptions.  
  - Sheet4–16 present derived intermediate datasets, key statistical outputs, and essential R code snippets used in the analyses. These sheets are intended to facilitate interpretation and reproducibility of the R scripts.

- **`moss_diversity.R`** contains the complete R code used to reproduce all statistical analyses and visualizations presented in the manuscript draft. The script uses **`env.csv`** and **`moss.csv`** as the primary input datasets.
