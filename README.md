# tufa_moss_diversity

## Repository Overview

This repository contains R scripts and processed data used for the manuscript:

**“Patterns and drivers of tufa moss diversity across the alpine petrifying springs landscape gradient on the Tibetan Plateau”** (in review)

The repository ensures full transparency and reproducibility of all statistical analyses and visualizations presented in the manuscript.

---

## Repository Structure

1. **`env.csv`**  
   Primary dataset used for all statistical analyses and visualizations. Additional details are provided in **`3_Supporting Table.xlsx`**, Sheet2 ("Variables data").

   The dataset includes:

   - Plot identifiers (`Sites`, `Plots`)
   - Geographic coordinates (`Lng`, `Lat`)
   - Habitat environmental variables (`Alt`, `Temp`, `RH`, `pH`, `Cond`, `OC`, `TN`, `TP`)
   - Moss morphological and physiological attributes (`MTC`, `MTN`, `MTP`, `Chl_a`, `Chl_b`, `Chl_a.b`, `Chl_a.b.1`, `MDA`, `TSS`, `POD`, `SOD`, `CAT`)
   - Diversity metrics (`Richness`, `Shannon`, `Simpson`, `NMDS1`)
   - Derived principal component variables (`PC1_diversity`, `PC1_geospatial`, `PC1_meteorological`, `PC1_Water`, `PC1_substrate`, `PC1_moss`, `PC1_photosynthesis`, `PC1_osmotic`, `PC1_enzyme`)

   The proportion of variance explained by the first principal component and correlations between principal components and original variables are reported in **`3_Supporting Table.xlsx`**, Sheet4 ("PCA result").

   Detailed variable definitions and abbreviations are available in **`3_Supporting Table.xlsx`**, Sheet3 ("Data description").

2. **`moss.csv`**  
   Moss community abundance dataset (cover values used as abundance proxies).  
   Detailed descriptions are available in **`3_Supporting Table.xlsx`**, Sheet1 ("Tufa moss data").

   This dataset is primarily used for:
   - Community composition analyses
   - Beta diversity partitioning
   - Ordination analyses (NMDS, RDA)
   - Distance-based analyses (Mantel tests)

3. **`3_Supporting Table.xlsx`**

   - Sheet1–3: Moss abundance data, variable data, and detailed data descriptions  
   - Sheet4–16: Derived intermediate datasets, key statistical results, and essential R code excerpts used in the analyses  

   These supplementary sheets facilitate transparency and reproducibility.

4. **`moss_diversity.R`**  
   Complete R script used to reproduce all statistical analyses and figures presented in the manuscript draft.

   The script includes:

   - Principal component analysis (PCA)
   - Linear regression analyses
   - One-way analysis of variance (ANOVA)
   - Multi-model inference
   - Non-metric multidimensional scaling (NMDS)
   - Indicator species analysis
   - Redundancy analysis (RDA)
   - Beta diversity partitioning
   - Mantel tests
   - Transformation-based RDA (tb-RDA)
   - Hierarchical partitioning
   - Variation partitioning
   - Structural equation modeling (SEM)

   Primary input datasets:
   - `env.csv`
   - `moss.csv`

---

## Reproducibility

All analyses were conducted in R.

Key R packages used in the analyses include:

**Community ecology & ordination**
- `vegan`
- `adespatial`
- `labdsv`
- `factoextra`

**Statistical modeling**
- `lme4`
- `MuMIn`
- `car`
- `performance`
- `effects`
- `agricolae`
- `emmeans`
- `piecewiseSEM`

**Multivariate partitioning**
- `rdacca.hp`

**Visualization**
- `ggplot2`
- `ggtern`
- `ggrepel`
- `ggprism`
- `ggupset`
- `cowplot`
- `networkD3`

**Data manipulation**
- `tidyverse`
- `reshape2`
- `readxl`
- `Rmisc`

**Spatial analysis**
- `geosphere`

Random permutation procedures (e.g., Mantel tests, hierarchical partitioning) may yield slightly varying p-values due to stochastic resampling. Reported statistics correspond to the specified permutation settings in the script.

---

## Citation

If you use the data, scripts, or conceptual framework from this repository, please cite:

Zhou et al. (in review). *Patterns and drivers of tufa moss diversity across the alpine petrifying springs landscape gradient on the Tibetan Plateau.*

This citation will be updated with journal information and DOI once available.

---

## Contact

For questions regarding the data, statistical analyses, or reproducibility of results, please contact the corresponding author via the email address provided in the R script (`moss_diversity.R`).

Academic collaboration and inquiries are welcome.
