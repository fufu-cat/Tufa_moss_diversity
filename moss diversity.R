# Title: Patterns and drivers of tufa moss diversity across alpine petrifying springs landscape gradient on the Tibetan Plateau
# Authors: Yuanyuan Zhou, Qiang Wei, Xiaodan Liang, Chengyi Li, Hu Chen, Zhihui Wang, Zhaohui Zhang*
# R code implementation and inspection: Yuanyuan Zhou (1181559954@qq.com); and Qiang Wei (wq_tendermoments@163.com)
# Date: February 2026
# Correspondence: Zhaohui Zhang (zhaozhang9@hotmail.com)




#### [REMARK ON VISUALIZATION] ####
# The spatial visualizations and maps presented in Figure 1 (a-d) were implemented using the ArcGIS 10.6 platform. This software was selected due to its superior precision in managing geographical coordinates and topographic data specific to tufa landscapes. Consequently, there is no corresponding R code for these map panels. 
# All other statistical plots (Figures 2-8) were initially generated as raw images using R 4.3.2.
# These raw outputs were subsequently refined, aestheticized, and assembled using Adobe Illustrator 2024 to ensure high-quality publication standards.
# Therefore, while the underlying data analysis and primary plotting are reproducible via this script, the final layout reflects post-processing optimization.
# Full names, abbreviations, and descriptions of all variables are provided in the Supporting Table (.xlsx): "Sheet3. Data description".
# data_desc <- read_excel("3_Supporting Table.xlsx", sheet = "3. Data description"); print(as_tibble(data_desc))




# Please ensure an active internet connection for the first run to install missing packages.
if (!require("pacman")) install.packages("pacman")
pacman::p_load( readxl, tidyverse,  Rmisc, vegan, reshape2, rdacca.hp, FactoMineR, labdsv, adespatial,
               factoextra, lme4, piecewiseSEM, MuMIn, car, performance, effects, ggtern, ggrepel,
               agricolae, emmeans, ggplot2, ggupset, cowplot, ggprism, linkET, networkD3, geosphere)




#### FIGURE 2: Effects of tufa landscapes on habitat environmental and moss ecophysiological traits ####
#### FIGURE 4a: Moss alpha diversity ####

moss_data <- read.csv("moss.csv")
library (vegan)
Richness <- rowSums(moss_data[,3:43] > 0)
Shannon <- diversity(moss_data[,3:43], "shannon", base = exp(1))
Simpson <- diversity(moss_data[,3:43], "simpson")
diversity <- data.frame(Richness, Shannon, Simpson)
# write.csv(diversity, file = "diversity.csv")

env <- read.csv("env.csv")
enviroment1 <- env %>% pivot_longer(cols = 6:27, names_to = "environment", values_to = "value")
enviroment1$environment <- factor(enviroment1$environment, levels = colnames(env)[6:27])
enviroment2 <- summarySE(enviroment1, measurevar = "value", groupvars = c("environment", "Sites"))
enviroment2$Sites <- factor(enviroment2$Sites, levels = c("A", "B", "C", "D", "E"))
# write.csv(enviroment2, "enviroment2.csv")

FIGURE_2_4a <- ggplot(data = enviroment2, aes(x = Sites, y = value, colour = Sites)) +
  geom_point(stat = "identity", size = 5, shape = 19) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se, x = Sites), linewidth = 0.7, width = 0.3, colour = "black") +
  geom_jitter(data = enviroment1, aes(colour = Sites), position = position_jitter(0), shape = 16, alpha = 0.2, size = 3) +
  facet_wrap(~ environment, scales = "free", ncol = 4) +
  theme_bw() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank(), axis.text.x = element_text(size = 10, colour = "black", angle = 0, hjust = 0.5), axis.text.y = element_text(size = 10)) +
  scale_colour_manual(values = c("#E3C847", "#FF9E6C", "#3399BC", "#9150B1", "#B83851")) +
  guides(colour = "none"); print(FIGURE_2_4a)
# ggsave(filename = "OLS.pdf", plot = FIGURE_2_4a, width = 19, height = 28, units = "cm", device = cairo_pdf)

names(env)
fit <- lm(Temp~ Sites, data = env) 
# fit <- lm(log(TN)~ Sites, data = env) # If residuals do not meet normality assumptions, consider log-transformation
summary(fit) # Model summary
check_normality(fit)
anova(fit) # Global significance test
LSD.test(fit, "Sites", console = TRUE) # Post-hoc grouping (a, b, c...)
# [Repeat the above steps manually for each variable from index 6 to 27]

## FIGURES 2 & 4a generation and statistical testing are now complete.




#### FIGURE 3a Relative abundance of dominant moss species ####

moss_bar <- read.csv("moss.csv")
moss_bar1 <- moss_bar[, c(1, 6, 16, 20, 25, 31, 34, 36)] # Seven dominant taxa were screened and selected for community composition analysis based on a relative abundance threshold of >1% and a relative frequency of >20%, as detailed in the Supporting Table (.xlsx): "1. Tufa moss data".
moss_bar3 <- reshape2::melt(moss_bar1, id.vars = 1, measure.vars = c(2:8), variable.name = "species", value.name = "value")
moss_bar4 <- summarySE(moss_bar3, measurevar = "value", groupvars = c("species", "Sites"))
# write.csv(moss_bar4, file = "moss_bar4.csv")

FIGURE_3a <- ggplot(data = moss_bar4, aes(x = Sites, y = value, fill = species)) +
  geom_bar(stat = "identity", position = "fill", color = "white", width = 0.8, linewidth = 0.5) +
  scale_fill_prism(palette = "summer") +
  coord_cartesian(ylim = c(0, 1), expand = TRUE, clip = 'off') +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.line = element_line(color = 'black'), axis.text = element_text(size = 12, color = "black"), axis.title = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 14, hjust = 0), legend.title = element_text(size = 10, color = "black"), legend.position = "right") +
  labs(x = "Sampling Sites", y = "Relative abundance (%)"); print(FIGURE_3a )
#ggsave('FIGURE_3a_Bar_Filtered.pdf',FIGURE_3a, width = 4, height = 3)

## FIGURE 3a generation is now complete.




#### FIGURE 3b Moss indicator species ####

moss_data <- read.csv("moss.csv")
set.seed(123)
indval_spe_obj <- indval(moss_data[, -(1:2)], moss_data$Sites, numitr = 10000)
indval_spe <- data.frame(
  maxcls    = indval_spe_obj$maxcls, # Index of the site with the highest indicator value
  indcls    = indval_spe_obj$indcls, # The maximum indicator value (IndVal) achieved across all sites
  pval      = indval_spe_obj$pval,   # Significance (p-value) derived from the permutation test
  frequency = apply(moss_data[, -(1:2)] > 0, 2, sum), # Total occurrences of each species
  relfrq    = indval_spe_obj$relfrq, # Matrix: Relative frequency of species in each site
  relabu    = indval_spe_obj$relabu, # Matrix: Relative abundance of species in each site
  indval    = indval_spe_obj$indval  # Matrix: Individual IndVal scores for all sites
)
indval_spe1 <- indval_spe %>% filter(pval < 0.05) %>% mutate(species = row.names(.),  maxcls = factor(maxcls, levels = c("1", "2", "3", "4", "5"), labels = c("A", "B", "C", "D", "E")))
# write.csv(indval_spe, "indval_spe.csv")

FIGURE_3b <- ggplot(indval_spe1, aes(x = maxcls, y = indcls, fill = maxcls, group = species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "white", linewidth = 0.5) +
  geom_text(aes(label = species, y = indcls/2), position = position_dodge(width = 0.9), vjust = 0.5, hjust = 0.5, angle = 90, size = 3, fontface = "italic", color = "white") +
  scale_fill_manual(values = c("#E3C847", "#FF9E6C", "#3399BC", "#9150B1", "#B83851")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11), axis.text = element_text(face = "bold", size = 11, color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Sampling Sites (Landscape Gradient)", y = "Indicator Value (IndVal)"); print(FIGURE_3b)
# ggsave("F_indval.pdf", FIGURE_3b, width = 12, height = 8, units = "cm")

## FIGURE 3b generation is now complete.




#### FIGURE 3c Shared and unique moss species ####

moss_upset <- read.csv("moss.csv")

moss_upset1 <- melt(moss_upset, id.vars = 1, measure.vars = 3:43, variable.name = "species", value.name = "value")
moss_upset1$value <- suppressWarnings(as.numeric(moss_upset1$value))
moss_upset2 <- summarySE(moss_upset1, measurevar = "value", groupvars  = c("species", "Sites"))
moss_upset3 <- moss_upset2$value > 0
moss_upset4 <- moss_upset2 %>% mutate(present = moss_upset3) %>% filter(present) %>% select(species, Sites)
moss_upset5 <- moss_upset4 %>% dplyr::group_by(species) %>% dplyr::summarise(Sites = list(Sites), .groups = "drop")
site_levels <- c("A", "B", "C", "D", "E")
moss_upset6 <- moss_upset5 %>% dplyr::rename(species = dplyr::any_of(c("species","Species","variable","name"))) %>% dplyr::mutate(Sites = lapply(Sites, function(x) site_levels[site_levels %in% as.character(x)]), sites_string = vapply(Sites, paste, character(1), collapse = ", ")) %>% dplyr::arrange(sites_string, species)
upset.csv <- moss_upset6 %>% dplyr::select(species, sites_string)
# write.csv(upset.csv, file = "upset.csv")

FIGURE_3c <- ggplot(moss_upset6, aes(x = Sites)) +
  scale_x_upset(sets = site_levels, order_by = "freq") +
  geom_bar(fill = "#88B766") +
  theme_combmatrix(combmatrix.panel.point.color.fill = "#3399BC", combmatrix.panel.point.color.empty = "grey90", combmatrix.panel.line.color = "#B83851", combmatrix.panel.line.size = 1.5, combmatrix.label.extra_spacing = 5, combmatrix.panel.point.size = 12) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = 1.2) +
  scale_y_continuous(breaks = seq(0, 9, by = 2), limits = c(0, 9), name = "Species richness (sp.)");  print(FIGURE_3c)
# ggsave('upset.pdf', FIGURE_3c, width = 6, height = 4) 

## FIGURE 3c generation is now complete.




#### FIGURE 4b-c Relative effects of environmental conditions and moss traits on community diversity ####

env <- read.csv("env.csv")
env_bio <- env[, c("Richness", "Shannon", "Simpson")]
cor(env_bio, method = "pearson")
env_bio_pca <- PCA(env_bio, scale.unit = TRUE, graph = FALSE)
summary(env_bio_pca)
get_eigenvalue(env_bio_pca)
pca_scores <- get_pca_ind(env_bio_pca)$coord
pcadata <- data.frame(PC1_diversity = pca_scores[, 1])
cor.test(env_bio$Richness, pca_scores[, 1], method = "pearson")
cor.test(env_bio$Shannon,  pca_scores[, 1], method = "pearson")
cor.test(env_bio$Simpson,  pca_scores[, 1], method = "pearson")
cor.test(env$NMDS1, pca_scores[, 1], method = "pearson") 
# Note: to integrate the three α-diversity indices (Richness, Shannon, and Simpson), we performed a principal component analysis (PCA) on standardized variables and extracted the first principal component (PC1) as a composite representation of moss α-diversity. PC1 was significantly correlated with each individual metric, indicating that it captures the shared variance among these indices. Accordingly, PC1 scores were used as the proxy for α-diversity in subsequent analyses (see Materials and Methods and Supporting Table Sheet4, “PCA result”).

env_scale <- as.data.frame(scale(env[, 3:37]))
options(na.action = "na.fail")
mul_diversity51 <- lm(PC1_diversity ~ Lng + Lat + Alt + Temp + RH + pH + Cond + OC + TN + TP, data = env_scale) # Full model
check_normality(mul_diversity51)
car::vif(mul_diversity51) 
car::vif(lm(PC1_diversity ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale)) # Repeat this line and remove the highest-VIF term each time (Lat, Alt, OC ...)
mul_diversity52 <- lm(PC1_diversity ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale) # Candidate model with VIF < 10
summary(mul_diversity52)
check_normality(mul_diversity52)
mul_selection52 <- MuMIn::dredge(mul_diversity52) # Model selection via model averaging (ΔAICc < 2) 
subset(mul_selection52, delta < 2)  # show best subset (2 models)
best52 <- MuMIn::model.avg(mul_selection52, subset = delta < 2)
summary(best52)  # averaged estimates & p-values
mul_diversity_best53 <- lm(PC1_diversity ~ Lng + Temp + RH, data = env_scale) # Refit a final model (for R2 reporting)
summary(mul_diversity_best53)
check_normality(mul_diversity_best53)
r2(mul_diversity_best53)  # for R2
AICc(mul_diversity51, mul_diversity52, mul_diversity_best53)

coef_tab <- as.data.frame(summary(best52)$coefficients)
coef_used <- coef_tab["full", , drop = FALSE]
aaa5 <- data.frame(variable = names(coef_used), Estimate = as.numeric(coef_used[1, ]), stringsAsFactors = FALSE)
ci_tab <- as.data.frame(confint(best52, full = TRUE))
ci_tab$variable <- rownames(ci_tab)
aaa5 <- aaa5 %>% dplyr::left_join(ci_tab, by = "variable") %>% dplyr::rename(lower = `2.5 %`, upper = `97.5 %`) %>% dplyr::filter(variable != "(Intercept)"); aaa5
aaa5$variable <- factor(aaa5$variable, levels = c("RH", "Temp", "Lng"))
aaa5 <- aaa5 %>% dplyr::mutate(type = c("Geospatial", "Climatic", "Climatic")); aaa5

f_diversity <- ggplot(aaa5) +
  geom_hline(yintercept = 0, linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(x = variable, y = Estimate, color = type, shape = type), size = 4) +
  geom_errorbar(aes(x = variable, ymin = lower, ymax = upper, color = type), linewidth = 0.8, width = 0.25) +
  scale_color_manual(values = c("Geospatial" = "#3399BC", "Climatic" = "#B83851")) +
  scale_shape_manual(values = c("Geospatial" = 15, "Climatic" = 16)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 1.5)) +
  ylab("Parameter estimates") +
  xlab(NULL); f_diversity

Estimate5  <- abs(aaa5$Estimate)
Geospatia5 <- sum(Estimate5[aaa5$type == "Geospatial"]) / sum(Estimate5) * 100
Climatic5  <- sum(Estimate5[aaa5$type == "Climatic"])  / sum(Estimate5) * 100
bb53 <- data.frame(value = c(Geospatia5, Climatic5), variable = factor(c("Geospatial", "Climatic"), levels = c("Geospatial", "Climatic"))); bb53

f_bar_diversity <- ggplot(bb53, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Geospatial" = "#3399BC", "Climatic" = "#B83851")) +
  theme_classic() +
  theme(legend.position = "none", axis.line = element_line(linetype = "solid"), axis.text.y = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_blank()) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_diversity

FIGURE_4b <- cowplot::plot_grid(f_bar_diversity, f_diversity, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_4b

## FIGURE 4b generation is now complete.


mul_diversity61 <- lm(PC1_diversity ~ MTC + MTN + MTP + Chl_a + Chl_b + Chl_a.b + MDA + TSS + POD + SOD + CAT, data = env_scale)
check_normality(mul_diversity61)
car::vif(mul_diversity61)
car::vif(lm(PC1_diversity ~ MTC + MTP + Chl_a + Chl_b + MDA + TSS + POD + SOD + CAT, data = env_scale))
mul_diversity62 <- lm(PC1_diversity ~ MTC + MTP + Chl_a + Chl_b + MDA + TSS + POD + SOD + CAT, data = env_scale)
summary(mul_diversity62)
check_normality(mul_diversity62)
mul_selection62 <- MuMIn::dredge(mul_diversity62)
subset(mul_selection62, delta < 2)  
best62 <- MuMIn::model.avg(mul_selection62, subset = delta < 2)
summary(best62)
mul_diversity_best63 <- lm(PC1_diversity ~ MTC + MTP + Chl_a + MDA + SOD + CAT, data = env_scale)
summary(mul_diversity_best63)
check_normality(mul_diversity_best63)
r2(mul_diversity_best63)
AICc(mul_diversity61, mul_diversity62, mul_diversity_best63)

aaa6 <- summary(best62)$coefficients[1,] %>% as.data.frame() %>% mutate(lower = confint(best62, full = TRUE)[, 1], upper = confint(best62, full = TRUE)[, 2]) %>% slice(-1)
names(aaa6)[1] <- "Estimate"
aaa6$variable <- row.names(aaa6)
aaa6$variable <- factor(aaa6$variable, levels = c("CAT", "SOD", "MDA", "Chl_a", "MTC", "MTP"))
aaa6 <- mutate(aaa6, type = c("Photosynthesis", "Nutrients", "enzyme", "Nutrients", "enzyme", "Osmotic"))

f_diversity1 <- ggplot(aaa6) +
  geom_hline(aes(yintercept = 0), linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0.25) +
  scale_color_manual(values = c("#863198", "#D8692F", "#0077D4", "#408E09")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 0.8)) +
  ylab("Parameter estimates"); f_diversity1

Estimate6 <- abs(aaa6[, 1])
Nutrients6       <- sum(Estimate6[c(2, 4)]) / sum(Estimate6) * 100
Photosynthesis6  <- sum(Estimate6[c(1)])    / sum(Estimate6) * 100
Osmotic6         <- sum(Estimate6[c(6)])    / sum(Estimate6) * 100
enzyme6          <- sum(Estimate6[c(5, 3)]) / sum(Estimate6) * 100
bb61 <- c(Nutrients6, Photosynthesis6, Osmotic6, enzyme6) %>% as.data.frame()
bb62 <- c("Nutrients", "Photosynthesis", "Osmotic", "enzyme") %>% as.data.frame()
bb63 <- cbind(bb61, bb62)
names(bb63) <- c("value", "variable")
bb63$variable <- factor(bb63$variable, levels = c("Nutrients", "Photosynthesis", "Osmotic", "enzyme"))

f_bar_diversity1 <- ggplot(bb63, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c(Nutrients = "#D8692F", Photosynthesis = "#408E09", Osmotic = "#0077D4", enzyme = "#863198")) +
  theme_classic() +
  theme(legend.position = "none", axis.line = element_line(linetype = "solid"), axis.text.y = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_text(size = 12, color = "black", face = "bold", hjust = 1), title = element_text(size = 12)) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_diversity1

FIGURE_4c <- cowplot::plot_grid(f_bar_diversity1, f_diversity1, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_4c

## FIGURE 4c generation is now complete.




#### FIGURE 5a-d Relative effects of tufa habitat conditions on moss traits ####

env <- read.csv("env.csv")
env_mos <- env[, c("MTC", "MTN", "MTP")] 
cor(env_mos, method = "pearson")  
env_mos_pca <- PCA(env_mos, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(env_mos_pca)
PC1_moss <- get_pca_ind(env_mos_pca)$coord[, 1]
pcadata <- cbind(pcadata, PC1_moss)
cor.test(env_mos$MTC, PC1_moss, method = "pearson", conf.level = 0.95)
cor.test(env_mos$MTN, PC1_moss, method = "pearson", conf.level = 0.95)
cor.test(env_mos$MTP, PC1_moss, method = "pearson", conf.level = 0.95)
# Mantel tests revealed highly significant correlations among MTC, MTN, and MTP, indicating strong collinearity. To reduce redundancy and avoid multicollinearity in subsequent models, we used the first principal component (PC1) derived from PCA of standardized variables as a composite stoichiometric proxy. Detailed procedures are described in the Materials and Methods section.

names(env_scale)
mul_nutrients11 <- lm(PC1_moss ~ Lng + Lat + Alt + Temp + RH + pH + Cond + OC + TN + TP, data = env_scale)
summary(mul_nutrients11)
check_normality(mul_nutrients11)
vif(mul_nutrients11)
vif(lm(PC1_moss ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale))
mul_nutrients12 <- lm( PC1_moss ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale)
summary(mul_nutrients12)
check_normality(mul_nutrients12)
mul_selection12 <- dredge(mul_nutrients12)
subset(mul_selection12, delta < 2)  
best12 <- model.avg(mul_selection12, subset = delta < 2)
summary(best12)
mul_nutrients_best13 <- lm(PC1_moss ~ Lng + Temp + RH + pH + Cond, data = env_scale)
summary(mul_nutrients_best13)
check_normality(mul_nutrients_best13)
r2(mul_nutrients_best13)
AICc(mul_nutrients11, mul_nutrients12, mul_nutrients_best13)

aaa1 <- summary(best12)$coefficients[1, ] %>% as.data.frame() %>% mutate(lower = confint(best12, full = TRUE)[, 1], upper = confint(best12, full = TRUE)[, 2]) %>%  slice(-1)  # remove intercept
names(aaa1)[1] <- "Estimate"
aaa1$variable <- row.names(aaa1)
aaa1$variable <- factor(aaa1$variable, levels = c("Cond", "pH", "RH", "Temp", "Lng"))
aaa1 <- mutate(aaa1, type = c("Geospatial", "Water", "Climatic", "Climatic", "Water"))

f_nutrients <- ggplot(aaa1) +
  geom_hline(aes(yintercept = 0), linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0.25) +
  scale_color_manual(values = c("#B83851", "#3399BC", "#E3C847")) +
  scale_shape_manual(values = c(16, 15, 17)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  ylab("Parameter estimates") +
  xlab(NULL); f_nutrients

Estimate1   <- abs(aaa1[, 1])
Geospatial1 <- sum(Estimate1[c(1)])   / sum(Estimate1) * 100
Climatic1   <- sum(Estimate1[c(3:4)]) / sum(Estimate1) * 100
Water1      <- sum(Estimate1[c(2, 5)]) / sum(Estimate1) * 100
bb11 <- c(Geospatial1, Climatic1, Water1) %>% as.data.frame()
bb12 <- c("Geospatial", "Climatic", "Water") %>% as.data.frame()
bb13 <- cbind(bb11, bb12)
names(bb13) <- c("value", "variable")
bb13$variable <- factor(bb13$variable, levels = c("Geospatial", "Climatic", "Water"))

f_bar_nutrients <- ggplot(bb13, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(Geospatial = "#3399BC", Climatic = "#B83851", Water = "#E3C847")) +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line = element_line(linetype = "solid"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black", face = "bold", hjust = 1),
    title = element_text(size = 12)
  ) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_nutrients

FIGURE_5a <- cowplot::plot_grid(f_bar_nutrients, f_nutrients, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_5a

## FIGURE 5a generation is now complete.


env_pho <- env[, c("Chl_a", "Chl_b", "Chl_a.b")] 
cor(env_pho, method = "pearson")
env_pho_pca <- PCA(env_pho, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(env_pho_pca)
PC1_photosynthesis <- get_pca_ind(env_pho_pca)$coord[, 1]
pcadata <- cbind(pcadata, PC1_photosynthesis)
cor.test(env_pho$Chl_a,   PC1_photosynthesis, method = "pearson", conf.level = 0.95)
cor.test(env_pho$Chl_b,   PC1_photosynthesis, method = "pearson", conf.level = 0.95)
cor.test(env_pho$Chl_a.b, PC1_photosynthesis, method = "pearson", conf.level = 0.95)

names(env_scale)
mul_Photosynthetic21 <- lm(PC1_photosynthesis ~ Lng + Lat + Alt + Temp + RH + pH + Cond + OC + TN + TP, data = env_scale)
check_normality(mul_Photosynthetic21) 
car::vif(mul_Photosynthetic21)
car::vif(lm(PC1_photosynthesis ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale))
mul_photosynthesis22 <- lm(PC1_photosynthesis ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale)
summary(mul_photosynthesis22)
check_normality(mul_photosynthesis22)  # should pass after VIF filtering
mul_selection22 <- MuMIn::dredge(mul_photosynthesis22)
subset(mul_selection22, delta < 2)
best22 <- MuMIn::model.avg(mul_selection22, subset = delta < 2)
summary(best22)
mul_photosynthesis_best23 <- lm(PC1_photosynthesis ~ Lng + Temp + RH + TN + Cond, data = env_scale)
summary(mul_photosynthesis_best23)
check_normality(mul_photosynthesis_best23)
r2(mul_photosynthesis_best23)
AICc(mul_Photosynthetic21, mul_photosynthesis22, mul_photosynthesis_best23)

aaa2 <- summary(best22)$coefficients[1, ] %>% as.data.frame() %>% mutate(lower = confint(best22, full = TRUE)[, 1], upper = confint(best22, full = TRUE)[, 2]) %>% slice(-1)  # remove intercept
names(aaa2)[1] <- "Estimate"
aaa2$variable <- row.names(aaa2)
aaa2$variable <- factor(aaa2$variable, levels = c("TN", "Cond", "RH", "Temp", "Lng"))
aaa2 <- mutate(aaa2, type = c("Water", "Geospatial", "Climatic", "Climatic", "Substrate")); aaa2

f_Photosynthetic <- ggplot(aaa2) +
  geom_hline(aes(yintercept = 0), linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0.25) +
  scale_color_manual(values = c("#B83851", "#3399BC", "#9150B1", "#E3C847")) +
  scale_shape_manual(values = c(16, 15, 18, 17)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.5, 1.5)) +
  ylab("Parameter estimates") +
  xlab(NULL); f_Photosynthetic

Estimate2 <- abs(aaa2[, 1])
Water2      <- sum(Estimate2[c(1)])   / sum(Estimate2) * 100
Geospatial2 <- sum(Estimate2[c(2)])   / sum(Estimate2) * 100
Climatic2   <- sum(Estimate2[c(3:4)]) / sum(Estimate2) * 100
Substrate2  <- sum(Estimate2[c(5)])   / sum(Estimate2) * 100
bb21 <- c(Geospatial2, Climatic2, Water2, Substrate2) %>% as.data.frame()
bb22 <- c("Geospatial", "Climatic", "Water", "Substrate") %>% as.data.frame()
bb23 <- cbind(bb21, bb22)
names(bb23) <- c("value", "variable")
bb23$variable <- factor(bb23$variable, levels = c("Geospatial", "Climatic", "Water", "Substrate"))

f_bar_Photosynthetic <- ggplot(bb23, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c(Geospatial = "#3399BC", Climatic = "#B83851", Water = "#E3C847", Substrate = "#9150B1")) +
  theme_classic() +
  theme(legend.position = "none", axis.line = element_line(linetype = "solid"), axis.text.y = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_text(size = 12, color = "black", face = "bold", hjust = 1), title = element_text(size = 12)) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_Photosynthetic

FIGURE_5b <- cowplot::plot_grid(f_bar_Photosynthetic, f_Photosynthetic, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_5b

## FIGURE 5b generation is now complete.


env_osm <- env[, c("MDA", "TSS")]
cor(env_osm, method = "pearson")
env_osm_pca <- PCA(env_osm, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(env_osm_pca)
PC1_osmotic <- get_pca_ind(env_osm_pca)$coord[, 1]
pcadata <- cbind(pcadata, PC1_osmotic)
cor.test(env_osm$MDA, PC1_osmotic, method = "pearson", conf.level = 0.95)
cor.test(env_osm$TSS, PC1_osmotic, method = "pearson", conf.level = 0.95)

names(env_scale)
mul_osmotic31 <- lm(PC1_osmotic ~ Lng + Lat + Alt + Temp + RH + pH + Cond + OC + TN + TP, data = env_scale)
check_normality(mul_osmotic31)
car::vif(mul_osmotic31)
car::vif(lm(PC1_osmotic ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale))
mul_osmotic32 <- lm(PC1_osmotic ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale)
summary(mul_osmotic32)
check_normality(mul_osmotic32)
mul_selection32 <- MuMIn::dredge(mul_osmotic32)
subset(mul_selection32, delta < 2)
best32 <- MuMIn::model.avg(mul_selection32, subset = delta < 2)
summary(best32)
mul_osmotic_best33 <- lm(PC1_osmotic ~ Lng + Temp + RH + TP + Cond, data = env_scale)
summary(mul_osmotic_best33)
check_normality(mul_osmotic_best33)
r2(mul_osmotic_best33)
AICc(mul_osmotic31, mul_osmotic32, mul_osmotic_best33)

aaa3 <- summary(best32)$coefficients[1, ] %>% as.data.frame() %>% mutate(lower = confint(best32, full = TRUE)[, 1], upper = confint(best32, full = TRUE)[, 2]) %>% slice(-1)
names(aaa3)[1] <- "Estimate"
aaa3$variable <- row.names(aaa3)
aaa3$variable <- factor(aaa3$variable, levels = c("TP", "Cond", "RH", "Temp", "Lng"))
aaa3 <- mutate(aaa3, type = c("Geospatial", "Climatic", "Climatic", "Substrate", "Water")); aaa3 

f_osmotic <- ggplot(aaa3) +
  geom_hline(aes(yintercept = 0), linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0.25) +
  scale_color_manual(values = c("#B83851", "#3399BC", "#9150B1", "#E3C847")) +
  scale_shape_manual(values = c(16, 15, 18, 17)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1.5)) +
  ylab("Parameter estimates") +
  xlab(NULL); f_osmotic

Estimate3   <- abs(aaa3[, 1])
Water3      <- sum(Estimate3[c(5)])   / sum(Estimate3) * 100
Geospatia3  <- sum(Estimate3[c(1)])   / sum(Estimate3) * 100
Climatic3   <- sum(Estimate3[c(2:3)]) / sum(Estimate3) * 100
Substrate3  <- sum(Estimate3[c(4)])   / sum(Estimate3) * 100
bb31 <- c(Geospatia3, Climatic3, Water3, Substrate3) %>% as.data.frame()
bb32 <- c("Geospatial", "Climatic", "Water", "Substrate") %>% as.data.frame()
bb33 <- cbind(bb31, bb32)
names(bb33) <- c("value", "variable")
bb33$variable <- factor(bb33$variable, levels = c("Geospatial", "Climatic", "Water", "Substrate"))

f_bar_osmotic <- ggplot(bb33, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  scale_fill_manual(values = c(Geospatial = "#3399BC", Climatic   = "#B83851", Water      = "#E3C847", Substrate  = "#9150B1")) +
  theme_classic() +
  theme(legend.position = "none", axis.line = element_line(linetype = "solid"), axis.text.y = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_text(size = 12, color = "black", face = "bold", hjust = 1), title = element_text(size = 12)) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_osmotic

FIGURE_5c <- cowplot::plot_grid(f_bar_osmotic, f_osmotic, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_5c

## FIGURE 5c generation is now complete.


env_enz <- env[, c("POD", "SOD", "CAT")]
cor(env_enz, method = "pearson")
env_enz_pca <- PCA(env_enz, scale.unit = TRUE, graph = FALSE)
get_eigenvalue(env_enz_pca)
PC1_enzyme <- get_pca_ind(env_enz_pca)$coord[, 1]
pcadata <- cbind(pcadata, PC1_enzyme)
cor.test(env_enz$POD, PC1_enzyme, method = "pearson", conf.level = 0.95)
cor.test(env_enz$SOD, PC1_enzyme, method = "pearson", conf.level = 0.95)
cor.test(env_enz$CAT, PC1_enzyme, method = "pearson", conf.level = 0.95)

names(env_scale)
mul_enzyme41 <- lm(PC1_enzyme ~ Lng + Lat + Alt + Temp + RH + pH + Cond + OC + TN + TP, data = env_scale)
check_normality(mul_enzyme41)
vif(mul_enzyme41)
vif(lm(PC1_enzyme ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale))
mul_enzyme42 <- lm(PC1_enzyme ~ Lng + Temp + RH + pH + Cond + TN + TP, data = env_scale)
summary(mul_enzyme42)
check_normality(mul_enzyme42)
mul_selection42 <- dredge(mul_enzyme42)
subset(mul_selection42, delta < 2)  
best42 <- model.avg(mul_selection42, subset = delta < 2)
summary(best42)  # averaged estimates & p-values
mul_enzyme_best43 <- lm(PC1_enzyme ~ Temp + RH + pH + Cond + TP + TN, data = env_scale)
summary(mul_enzyme_best43)
check_normality(mul_enzyme_best43)
r2(mul_enzyme_best43)
AICc(mul_enzyme41, mul_enzyme42, mul_enzyme_best43)

aaa4 <- summary(best42)$coefficients[1, ] %>% as.data.frame() %>% mutate(lower = confint(best42, full = TRUE)[, 1], upper = confint(best42, full = TRUE)[, 2]) %>% slice(-1) 
names(aaa4)[1] <- "Estimate"
aaa4$variable <- row.names(aaa4)
aaa4$variable <- factor(aaa4$variable, levels = c("TN", "TP", "Cond", "pH", "RH", "Temp"))
aaa4 <- mutate(aaa4, type = c("Climatic", "Substrate", "Water", "Water", "Climatic", "Substrate")); aaa4

f_enzyme <- ggplot(aaa4) +
  geom_hline(aes(yintercept = 0), linewidth = 0.8, colour = "grey70", linetype = "dashed") +
  geom_point(aes(y = Estimate, x = variable, color = type, shape = type), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, x = variable, color = type), size = 0.8, width = 0.25) +
  scale_color_manual(values = c("#B83851", "#9150B1", "#E3C847")) +
  scale_shape_manual(values = c(16, 18, 17)) +
  coord_flip() +
  theme_classic() +
  theme(axis.line = element_line(linetype = "solid"), panel.background = element_rect(fill = NA), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12, face = "bold"), strip.text.x = element_text(size = 12), legend.position = "none") +
  scale_x_discrete(position = "top") +
  scale_y_continuous(expand = c(0, 0), limits = c(-0.8, 1.2)) +
  ylab("Parameter estimates") +
  xlab(NULL); f_enzyme

Estimate4  <- abs(aaa4[, 1])
Water4     <- sum(Estimate4[c(3:4)]) / sum(Estimate4) * 100
Climatic4  <- sum(Estimate4[c(1, 5)]) / sum(Estimate4) * 100
Substrate4 <- sum(Estimate4[c(2, 6)]) / sum(Estimate4) * 100
bb41 <- c(Climatic4, Water4, Substrate4) %>% as.data.frame()
bb42 <- c("Climatic", "Water", "Substrate") %>% as.data.frame()
bb43 <- cbind(bb41, bb42)
names(bb43) <- c("value", "variable")
bb43$variable <- factor(bb43$variable, levels = c("Climatic", "Water", "Substrate"))

f_bar_enzyme <- ggplot(bb43, aes(x = "", y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c(Climatic = "#B83851", Water = "#E3C847", Substrate = "#9150B1")) +
  geom_text(aes(label = scales::percent(value / sum(value), accuracy = 0.1)), position = position_fill(vjust = 0.5), colour = "white", size = 3, fontface = "bold") +
  theme_classic() +
  theme(legend.position = "none", axis.line = element_line(linetype = "solid"), axis.text.y = element_text(size = 12, color = "black", face = "bold"), axis.text.x = element_text(size = 12, color = "black", face = "bold", hjust = 1), title = element_text(size = 12)) +
  labs(x = "", y = "Relative effect of estimates (%)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)); f_bar_enzyme

FIGURE_5d <- cowplot::plot_grid(f_bar_enzyme, f_enzyme, ncol = 2, nrow = 1, rel_widths = c(1, 2.5)); FIGURE_5d

## FIGURE 5c generation is now complete.

FIGURE_5 <- plot_grid(FIGURE_5a, FIGURE_5b, FIGURE_5c, FIGURE_5d, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 14, label_x = 0.02, label_y = 0.98); FIGURE_5

## FIGURE 5 generation is now complete.




#### FIGURE 6 Structural equation model ####

env <- read.csv("env.csv", check.names = FALSE)
psem_data <- data.frame(env[, c(6:9, 29:30, 33:37)]) # due to different units among the first variables, we use scaled predictors in regressions
names(psem_data)
vif(lm(PC1_diversity ~ scale(Temp) + scale(RH) + scale(pH) + scale(Cond), data = psem_data)) # VIF check (all VIF < 3)
model_environment <- lm(PC1_diversity ~ scale(Temp) + scale(RH) + scale(pH) + scale(Cond), data = psem_data)
summary(model_environment)
check_normality(model_environment)
coefs(model_environment, standardize = "scale")  # standardized coefficients
b_Temp <- summary(model_environment)$coefficients["scale(Temp)", "Estimate"]; b_Temp
b_RH   <- summary(model_environment)$coefficients["scale(RH)",   "Estimate"]; b_RH
b_pH   <- summary(model_environment)$coefficients["scale(pH)",   "Estimate"]; b_pH
b_Cond <- summary(model_environment)$coefficients["scale(Cond)", "Estimate"]; b_Cond
new_environment <- b_Temp * scale(psem_data$Temp) + b_RH * scale(psem_data$RH) + b_pH * scale(psem_data$pH) + b_Cond * scale(psem_data$Cond) # Weighted composite (on scaled predictors)
psem_data$new_environment <- as.numeric(new_environment)
summary(lm(PC1_diversity ~ PC1_geospatial, data = psem_data))
coefs(lm(PC1_diversity ~ new_environment, data = psem_data))
# Note: Several predictors showed strong intercorrelations (see Mantel test results). To reduce redundancy and avoid multicollinearity while keeping the SEM parsimonious, we replaced highly correlated variables with principal-component proxies where appropriate and further constructed composite variables. Specifically, composites were calculated as weighted sums of standardized predictors, with weights derived from regression coefficients estimated in preliminary linear models.

psem_data1 <- data.frame(psem_data[, 5:12])
names(psem_data1)
mosses.list <- list(
  lm(PC1_substrate ~ PC1_geospatial + new_environment, data = psem_data1),
  lm(new_environment ~ PC1_geospatial, data = psem_data1),
  lm(PC1_photosynthesis ~ PC1_substrate + new_environment, data = psem_data1),
  lm(PC1_osmotic ~ PC1_substrate + new_environment, data = psem_data1),
  lm(PC1_moss ~ PC1_substrate + new_environment, data = psem_data1),
  lm(PC1_enzyme ~ PC1_substrate + new_environment, data = psem_data1),
  lm(PC1_diversity ~ PC1_geospatial + PC1_substrate + new_environment + PC1_moss + PC1_photosynthesis + PC1_osmotic + PC1_enzyme, data = psem_data1))
mosses.psem <- as.psem(mosses.list)
summary(mosses.psem)

mosses.psem1 <- update(mosses.psem, PC1_osmotic ~ PC1_geospatial + PC1_substrate + new_environment) # Stepwise updates (add missing paths suggested by diagnostics)
summary(mosses.psem1)
mosses.psem2 <- update(mosses.psem1, PC1_enzyme ~ PC1_geospatial + PC1_substrate + new_environment)
summary(mosses.psem2)
mosses.psem3 <- update(mosses.psem2, PC1_enzyme ~ PC1_moss + PC1_geospatial + PC1_substrate + new_environment)
summary(mosses.psem3)
mosses.psem4 <- update(mosses.psem3, PC1_photosynthesis ~ PC1_geospatial + PC1_substrate + new_environment)
summary(mosses.psem4)

plot(mosses.psem4)
sem_coefs <- coefs(mosses.psem4, standardize = "scale"); sem_coefs
bbb <- rsquared(mosses.psem4); bbb
AIC(mosses.psem4)
fisherC(mosses.psem4)

# All statistical results underlying Figure 6 were obtained from the final piecewise structural equation model described above. The graphical representation of the final SEM was produced and refined in Adobe Illustrator 2024 to ensure clarity and visual consistency.




#### FIGURE 7a Moss beta diversity ####

NMDS_data <- read.csv("moss.csv", check.names = FALSE)
names(NMDS_data)
set.seed(1234)
NMDS <- metaMDS(NMDS_data[, 3:43], distance = "bray", k = 2, trymax = 200)
NMDS$stress
stressplot(NMDS)
NMDS_points <- as.data.frame(NMDS$points)
colnames(NMDS_points)[1:2] <- c("NMDS1", "NMDS2")
NMDS_result <- NMDS_data %>% dplyr::select(Sites, Plots) %>% dplyr::bind_cols(NMDS_points); NMDS_result
# write.csv(NMDS_result, "NMDS.csv", row.names = FALSE)
centroids <- NMDS_result %>% dplyr::group_by(Sites) %>% dplyr::summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop"); centroids

FIGURE_7a <- ggplot(NMDS_result, aes(NMDS1, NMDS2, color = Sites)) +
  geom_point(size = 4) +
  scale_color_manual(values = c(A = "#E3C847", B = "#FF9E6C", C = "#3399BC", D = "#9150B1", E = "#B83851")) +
  geom_segment( data = NMDS_result %>% left_join(centroids, by = "Sites", suffix = c("", "_centroid")), aes(x = NMDS1_centroid, y = NMDS2_centroid, xend = NMDS1, yend = NMDS2), size = 0.6, alpha = 0.3) +
  geom_point(data = centroids, aes(x = NMDS1, y = NMDS2, color = Sites), shape = 3, size = 3, stroke = 1) +
  stat_ellipse(geom = "polygon", level = 0.5, linetype = 2, size = 0.6, aes(color = Sites), fill = NA, alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "solid", size = 0.6, color = "grey", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid", size = 0.6, color = "grey", alpha = 0.5) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11), axis.text  = element_text(face = "bold", size = 11, color = "black", hjust = 0.5), panel.grid = element_blank(), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none"); FIGURE_7a

# ggsave("NMDS.pdf", FIGURE_7a, width = 9, height = 8, units = "cm")

## FIGURE 7a generation is now complete.




#### FIGURE 7b Ternary plot of moss beta diversity partitioning ####

abund <- read.csv("moss.csv", check.names = FALSE)
beta_jac <- beta.div.comp(abund[, 3:43], coef = "J", quant = FALSE)
dist_jac_D <- beta_jac$D
dist_jac_repl <- beta_jac$repl
dist_jac_rich <- beta_jac$rich
jac_D    <- as.vector(dist_jac_D)
jac_repl <- as.vector(dist_jac_repl)
jac_rich <- as.vector(dist_jac_rich)
tern_df <- data.frame(Similarity = 1 - jac_D, Turnover   = jac_repl, Nestedness = jac_rich); tern_df
centroid_df <- data.frame(Similarity = mean(tern_df$Similarity, na.rm = TRUE), Turnover   = mean(tern_df$Turnover,   na.rm = TRUE), Nestedness = mean(tern_df$Nestedness, na.rm = TRUE)); centroid_df

FIGURE_7b <- ggtern(tern_df, aes(x = Similarity, y = Turnover, z = Nestedness)) +
  geom_point(alpha = 0.55, size = 3.2, shape = 16, colour = "#B83851") +
  geom_point(data = centroid_df, aes(x = Similarity, y = Turnover, z = Nestedness),
             colour = "black", size = 4.5, shape = 16) +
  labs(x = "Similarity (1 - βjac)", y = "Turnover (βrepl)", z = "Nestedness (βrich)") +
  theme_rgbw() +
  theme(legend.position = "none"); FIGURE_7b

# ggsave("beta_ternary.pdf", FIGURE_7b, width = 10, height = 10, units = "cm")

## FIGURE 7b generation is now complete.



#### FIGURE 7c-d Moss community distance relationships ####

env   <- read.csv("env.csv",  check.names = FALSE)
abund <- read.csv("moss.csv", check.names = FALSE)
dist_bray <- vegdist(abund[, 3:43], method = "bray")
env_vars <- c("Temp", "Alt", "pH")
dist_env <- dist(scale(env[, env_vars]), method = "euclidean")
eco_vars <- c("MTC", "MTN", "Chl_a", "Chl_b", "MDA", "POD", "SOD")
dist_eco <- dist(scale(env[, eco_vars]), method = "euclidean")
coords <- env[, c("Lng", "Lat")]
dist_geo <- as.dist(distm(coords, fun = distHaversine))

dist_vec <- function(d, div = 1) as.vector(d) / div
reg_plot <- function(df, x, y, xlab, ylab, point_col = "#3399BC") {
  ggplot(df, aes_string(x = x, y = y)) +
    geom_point(size = 3.2, alpha = 0.6, colour = point_col, shape = 16) +
    geom_smooth(method = "lm", colour = "black", alpha = 0.25, linewidth = 0.9) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme( panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"), axis.text.x = element_text(face = "bold", colour = "black", size = 11), axis.text.y = element_text(face = "bold", colour = "black", size = 11), axis.title  = element_text(face = "bold", colour = "black", size = 13))
}

df_dist <- data.frame(geo_km = dist_vec(dist_geo, div = 1000), env = dist_vec(dist_env), eco = dist_vec(dist_eco), bray = dist_vec(dist_bray)); df_dist
p_geo_bray <- reg_plot(df_dist, "geo_km", "bray", "Geographic distance (km)", "Bray-Curtis dissimilarity", point_col = "#9150B1"); p_geo_bray
a1 <- mantel(dist_bray, dist_geo, method = "spearman", permutations = 9999, na.rm = TRUE); a1 # Because the Mantel test is based on permutation procedures, p-values may vary slightly among runs due to stochastic resampling.
p_env_bray <- reg_plot(df_dist, "env", "bray", "Environmental distance", "Bray-Curtis dissimilarity", point_col = "#3399BC"); p_env_bray
a2 <- mantel(dist_bray, dist_env, method = "spearman", permutations = 9999, na.rm = TRUE); a2
p_eco_bray <- reg_plot(df_dist, "eco", "bray", "Trait distance", "Bray-Curtis dissimilarity", point_col = "#B83851"); p_eco_bray
a3 <- mantel(dist_bray, dist_eco, method = "spearman", permutations = 9999, na.rm = TRUE); a3
p_geo_eco  <- reg_plot(df_dist, "geo_km", "eco", "Geographic distance (km)", "Trait distance", point_col = "#9150B1"); p_geo_eco
a4 <- mantel(dist_geo, dist_eco, method = "spearman", permutations = 9999, na.rm = TRUE); a4
p_env_eco  <- reg_plot(df_dist, "env", "eco", "Environmental distance", "Trait distance", point_col = "#3399BC"); p_env_eco
a5 <- mantel(dist_env, dist_eco, method = "spearman", permutations = 9999, na.rm = TRUE); a5

FIGURE_7cd <- plot_grid(p_geo_bray, p_env_bray, p_eco_bray, p_geo_eco, p_env_eco, ncol = 3, nrow = 2, align = "hv"); FIGURE_7cd 

# ggsave("distance_line.pdf", FIGURE_7cd, width = 21, height = 13.5, units = "cm")

## FIGURE 7c-d generation is now complete.




#### FIGURE 8a Mantel test analysis of drivers of community beta diversity ####

env  <- read.csv("env.csv")
moss <- read.csv("moss.csv")
moss_mantel <- moss[, -c(1, 2)]  
ST  <- moss_mantel[, c(1:3, 6, 11, 13, 15:21)]
TT  <- moss_mantel[, c(4:5, 7:10, 12, 14, 22:23)]
W   <- moss_mantel[, c(24:41)]
Dom <- moss_mantel[, c(3:4, 12, 14, 17:18, 20, 23, 28:29, 32, 34, 38, 39)]
moss_mantel_data <- data.frame(ST, TT, W, Dom)
names(moss_mantel_data)  
any(duplicated(names(moss_mantel_data))) 
bin_mantel <- function(df) {
  df %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.3, Inf), labels = c("< 0.2", "0.2-0.3", ">= 0.3")), pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05")))
}
spec_map <- list(ST = 1:13, TT = 14:23, W = 24:41, Dom = 42:55)
mantel_env1 <- mantel_test(moss_mantel_data, env[, 3:12], spec_select = spec_map) %>% bin_mantel()
# write.csv(mantel_env1, file = "mantel.csv", row.names = FALSE)

F_mantel1 <- qcorrplot(correlate(env[, 3:12], method = "pearson"), type = "lower", diag = FALSE) +
  geom_square() +
  geom_mark(sep = "\n", sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, color = "white", size = 3) +
  geom_diag_label(mapping = aes(x = .x + 1), hjust = 2) +
  geom_couple(data = mantel_env1, aes(xend = .xend + 1, colour = pd, size = rd), alpha = 0.9, node.colour = c("red", "grey80"), node.fill   = c("red", "grey80"), node.size   = c(7, 4)) +
  scale_fill_gradientn(colors = c("#CE4368", "white", "#548235"), limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
  scale_size_manual(values = c("< 0.2" = 0.1, "0.2-0.3" = 1.5, ">= 0.3" = 3)) +
  scale_colour_manual(values = c("< 0.001" = "#E39C58", "0.001 - 0.01" = "#88B766", "0.01 - 0.05" = "#9999CC", ">= 0.05" = "grey80")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey45"), order = 2), colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1), fill   = guide_colorbar(title = "Pearson's r", order = 3)); F_mantel1

# ggsave("mental.pdf", F_mantel1, width = 8.5, height = 8)


mantel_env2 <- mantel_test(moss_mantel_data, env[, 13:24], spec_select = spec_map) %>% bin_mantel()
# write.csv(mantel_env2, file = "mantel1.csv", row.names = FALSE)

F_mantel2 <- qcorrplot(correlate(env[, 13:24], method = "pearson"), type = "upper", diag = FALSE) +
  geom_square() +
  geom_mark(sep = "\n", sig_level = c(0.05, 0.01, 0.001), sig_thres = 0.05, color = "white", size = 3) +
  geom_diag_label(mapping = aes(x = .x + 1), hjust = 1) +
  geom_couple(data = mantel_env2, aes(xend = .xend - 1.5, colour = pd, size = rd), alpha = 0.9, node.colour = c("red", "grey80"), node.fill   = c("red", "grey80"), node.size   = c(7, 4)) +
  scale_fill_gradientn(colors = c("#551461", "#F6F6F6", "#01451B"), limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
  scale_size_manual(values = c("< 0.2" = 0.1, "0.2-0.3" = 1.5, ">= 0.3" = 3)) +
  scale_colour_manual(values = c( "< 0.001" = "#E39C58", "0.001 - 0.01" = "#88B766", "0.01 - 0.05" = "#9999CC", ">= 0.05" = "grey80")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey45"), order = 2), colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1), fill   = guide_colorbar(title = "Pearson's r", order = 3)); F_mantel2

# ggsave("mental2.pdf", F_mantel2, width = 8.5, height = 8)
# Note: Adobe Illustrator 2024 was used for final figure polishing and layout.

## FIGURE 8a generation is now complete.




#### FIGURE 8b RDA and hierarchical partitioning reveal key drivers of moss communities ####

env <- read.csv("env.csv")
moss_RDA <- read.csv("moss.csv")
moss_RDA <- moss_RDA[, -c(1, 2)]
# Note: Many species contain excessive zeros; they cluster around the origin with weak environmental preference. Here we use species with relative abundance > 1% (and you may alternatively test stricter filters; results are similar).
moss_RDA1 <- moss_RDA[, c(3:4, 12, 14, 17:18, 20, 23, 28:30, 32, 34, 38:39)]
# moss_RDA1 <- moss_RDA[, c(4, 14, 18, 23, 29, 32, 34)] # alternative (stricter) set
# moss_RDA1 <- moss_RDA # all species (optional)
env_RDA1 <- env[, 30:37] # Explanatory variables (already PC1 proxies; no unit differences)
names(moss_RDA1)
names(env_RDA1)
moss_RDA2 <- decostand(moss_RDA1, method = "hellinger")
dca <- decorana(moss_RDA2); dca 
tb_RDA <- rda(moss_RDA2 ~ ., data = env_RDA1, scale = FALSE)
vif.cca(tb_RDA)  # VIF < 10 indicates low multicollinearity
tb_RDA_sum <- summary(tb_RDA, scaling = 1); tb_RDA_sum
tb_RDA_test_global <- anova.cca(tb_RDA, permutations = 999); tb_RDA_test_global # global test
tb_RDA_test_axis <- anova.cca(tb_RDA, by = "axis", permutations = 999); tb_RDA_test_axis # axis-by-axis
tb_RDA_envfit <- envfit(tb_RDA, env_RDA1, permutations = 999); tb_RDA_envfit

tb_RDA_site <- data.frame(tb_RDA_sum$site[, 1:2])
tb_RDA_site$Sites <- env$Sites
tb_RDA_env <- data.frame(tb_RDA_sum$biplot[, 1:2])
tb_RDA_env$env <- rownames(tb_RDA_env)
tb_RDA_species <- data.frame(tb_RDA_sum$species[, 1:2])
tb_RDA_species$species <- rownames(tb_RDA_species)

f_RDA <- ggplot(tb_RDA_site, aes(RDA1, RDA2)) +
  geom_point(aes(color = Sites), size = 5) +
  scale_colour_manual(values = c(A = "#E3C847", B = "#FF9E6C", C = "#3399BC", D = "#9150B1", E = "#B83851")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_segment( data = tb_RDA_env, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 22), colour = "black", linewidth = 0.5) +
  geom_text_repel(data = tb_RDA_env, aes(x = RDA1, y = RDA2, label = env), size = 4, segment.colour = "black") +
  geom_segment(data = tb_RDA_species, aes(x = 0, y = 0, xend = RDA1 * 0.5, yend = RDA2 * 0.5), arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 22), colour = "red", linewidth = 0.5) +
  geom_text_repel(data = tb_RDA_species, aes(x = RDA1, y = RDA2, label = species), size = 4, segment.colour = "black") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11), axis.text  = element_text(face = "bold", size = 11, color = "black", hjust = 0.5), panel.grid = element_blank(), plot.background = element_blank(), legend.position = "none") +
  xlab("RDA1 (26.35%)") + ylab("RDA2 (12.95%)"); f_RDA

# ggsave("RDA.pdf", f_RDA, width = 9, height = 8)


RDA_hp <- rdacca.hp(moss_RDA2, env_RDA1, method = "RDA", type = "R2", scale = FALSE) # Individual effects (R2-based)
set.seed(123)
permu_hp <- permu.hp(dv = moss_RDA2, iv = env_RDA1, method = "RDA", type = "R2", permutations = 999) 
hp_result <- data.frame(RDA_hp$Hier.part, p_value = permu_hp[, 2])
hp_result$env <- rownames(hp_result)
hp_result$env <- factor(hp_result$env, levels = c("PC1_substrate", "PC1_photosynthesis", "PC1_Water", "PC1_meteorological", "PC1_moss", "PC1_osmotic", "PC1_enzyme", "PC1_geospatial")); hp_result
# Because significance was assessed via permutation, p-values may vary slightly across runs. The estimated unique (individual) contributions are deterministic given the same input data and model settings.

f_HP <- ggplot(hp_result, aes(x = env, y = Individual, fill = env)) +
  geom_col(width = 0.6, color = "gray") +
  labs(x = "", y = "Individual effect") +
  scale_fill_prism(palette = "summer") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", color = "black", size = 11), axis.text  = element_text(face = "bold", size = 11, color = "black", hjust = 0.5), panel.grid = element_blank(), plot.background = element_blank(), legend.position = "none") +
  coord_flip(); f_HP

#ggsave("hierarchical partitioning.pdf", f_HP, width = 7, height = 8)

FIGURE_8b <- plot_grid(f_RDA, f_HP, ncol = 2); FIGURE_8b

## FIGURE 8b generation is now complete.




#### FIGURE 8c Venn diagram from variance partitioning analysis ####

envir <- list( Biotic= env_RDA1[c("PC1_moss","PC1_osmotic" ,"PC1_enzyme" )],env333=env_RDA1[c("PC1_Water","PC1_meteorological")],Abiotic=env_RDA1[c("PC1_geospatial")])
RDA_hp1 <- rdacca.hp(moss_RDA2, envir, method = 'RDA', type = "R2", scale = FALSE); RDA_hp1
permu_hp <- permu.hp(dv = moss_RDA2, iv = envir, method = 'RDA', type = 'R2', permutations = 999)
rda.vpa <- varpart(moss_RDA2, ~PC1_moss + PC1_osmotic + PC1_enzyme, ~PC1_Water+PC1_meteorological,~PC1_geospatial,data = env_RDA1,transfo = "hel", chisquare = FALSE)
dssa <- rda.vpa$part$contr1; rda.vpa
# write.csv(dssa,"dssa.csv")
anova.cca(rda(moss_RDA2~PC1_Water + PC1_meteorological+Condition(PC1_moss + PC1_osmotic + PC1_geospatial+ PC1_enzyme),env_RDA1),permutations=999)
anova.cca(rda(moss_RDA2~PC1_geospatial+Condition(PC1_moss + PC1_osmotic + PC1_Water + PC1_meteorological+ PC1_enzyme),env_RDA1),permutations=999)
anova.cca(rda(moss_RDA2~PC1_moss + PC1_osmotic +PC1_enzyme+Condition( PC1_Water + PC1_meteorological+ PC1_geospatial),env_RDA1),permutations=999)

showvarparts(3, bg=2:4)
plot(rda.vpa, cutoff = -Inf, cex = 1, digits = 2, bg = 2:5, Xnames=c("eco","env","geo")) 

## FIGURE 8c generation is now complete.


#### All statistical analyses and visualizations presented in the main text have been completed ####
#END



