# =========================================
# Biochar Cocoa Farm Analysis Script
# =========================================
# setwd("/Users/otissenalor/Documents/Otis Masters Thesis") # If you have a project, you don't need to set absolute paths.
# THIS AVOIDS PROBLEMS IN THE SECOND LINE OF CODE ALREADY. THAT IS THE BEAUTY OF THE Rproj. ;)
# -----------------------------
# 1. Setup
# -----------------------------
if(!require("Require")) install.packages("Require")
Require::Require(c("MASS","ggplot2","dplyr","lme4","car","emmeans","multcomp","tidyr", "ggpubr", "glmmTMB"))

# Consistent ggplot theme
theme_set(theme_minimal(base_size = 14))

# -----------------------------
# 2. Load Data
# -----------------------------
df <- read.csv("Updatedcsv2.csv", stringsAsFactors = TRUE)

# Create unique TreeID
df <- df %>%
  mutate(TreeUIDs = interaction(Farm, TreeUID, drop = TRUE))

# Quick check
str(df)
summary(df)

# -----------------------------
# 3. NDVI Analysis (Farm-Level Δ)
# -----------------------------
ndvi_df <- df %>%
  group_by(Farm, Treatment, Season) %>%
  summarise(NDVI = mean(Mean.NDVI, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Season, values_from = NDVI) %>%
  mutate(Delta = Dry - Wet)   # Δ = Dry − Wet
# --> This is very risky! If your data is not correctly organized, you will get the wrong results
# The correct way to do this is having a unique ID for each tree (i.e., each repeated measure )

# Boxplot of Δ NDVI by Treatment
ggplot(ndvi_df, aes(x = Treatment, y = Delta, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "Δ NDVI (Dry - Wet) by Treatment", y = "Δ NDVI") +
  theme(legend.position = "none")
# --> I love the plot but would be nice to see the significancy

# Clean data: remove NA and ensure Treatment is a factor
ndvi_df_clean <- ndvi_df %>%
  filter(!is.na(Delta), !is.na(Treatment)) %>%
  mutate(Treatment = factor(Treatment))

# Define pairwise comparisons
ndvi_df_clean2 <- ndvi_df_clean %>% filter(Treatment != "TS_O")

my_comparisons <- list(
  c("I_BC", "I_BCM"),
  c("I_BC", "I_BCO"),
  c("TS_BC", "TS_BCM"),
  c("TS_BC", "TS_BCO")
)


# Boxplot with ANOVA and pairwise significance
ggplot(ndvi_df_clean, aes(x = Treatment, y = Delta, fill = Treatment)) +
  geom_boxplot() +
  labs(title = "Δ NDVI (Dry - Wet) by Treatment", y = "Δ NDVI") +
  theme(legend.position = "none") +
  # ANOVA p-value at top
  stat_compare_means(method = "anova",
                     label.y = max(ndvi_df_clean$Delta, na.rm = TRUE) + 0.05) +
  # Pairwise t-tests
  stat_compare_means(comparisons = my_comparisons,
                     method = "t.test",
                     label = "p.signif",
                     label.y = max(ndvi_df_clean$Delta, na.rm = TRUE) +
                       seq(0.1, 0.1*length(my_comparisons), by = 0.05))


# Exclude Farm 1 from the dataset
ndvi_df_sub <- subset(ndvi_df, Farm != "Farm1")

# Run contrasts on the filtered data
contrast_TSBCM_TSM <- run_contrast(ndvi_df_sub, "TS_BCM", "TS_M")
contrast_IBCM_TSBCM <- run_contrast(ndvi_df_sub, "I_BCM", "TS_BCM")
contrast_TSBCO_TSO <- run_contrast(ndvi_df_sub, "TS_BCO", "TS_O")
contrast_IBCO_TSBCO <- run_contrast(ndvi_df_sub, "I_BCO", "TS_BCO")

# -----------------------------
# 4. NDVI Contrasts & Paired t-tests
# -----------------------------
run_contrast <- function(df, trt1, trt2){
  contrast <- df %>%
    filter(Treatment %in% c(trt1, trt2)) %>%
    dplyr::select(Farm, Treatment, Delta) %>%
    tidyr::pivot_wider(names_from = Treatment, values_from = Delta) %>%
    mutate(Diff = !!sym(trt1) - !!sym(trt2))

  t_res <- t.test(contrast[[trt1]], contrast[[trt2]], paired = TRUE)
  list(contrast = contrast, t_test = t_res)
}

# ---->  Did you check if your data allows you to perform parametric stats?


# Run contrasts
contrast_TSBCM_TSM <- run_contrast(ndvi_df, "TS_BCM", "TS_M")
contrast_IBCM_TSBCM <- run_contrast(ndvi_df, "I_BCM", "TS_BCM")

# Exclude Farm 1 before running the contrast
ndvi_df_sub <- subset(ndvi_df, Farm != "Farm1")

# Now run the contrast only on farms where both TS_BCO and TS_O exist
contrast_TSBCO_TSO <- run_contrast(ndvi_df_sub, "TS_BCO", "TS_O")

contrast_IBCO_TSBCO <- run_contrast(ndvi_df, "I_BCO", "TS_BCO")

# Summarize NDVI contrast results
ndvi_results <- data.frame(
  Contrast = c("TS_BCM vs TS_M", "I_BCM vs TS_BCM", "TS_BCO vs TS_O", "I_BCO vs TS_BCO"),
  Mean_Diff = c(mean(contrast_TSBCM_TSM$contrast$Diff, na.rm = TRUE),
                mean(contrast_IBCM_TSBCM$contrast$Diff, na.rm = TRUE),
                mean(contrast_TSBCO_TSO$contrast$Diff, na.rm = TRUE),
                mean(contrast_IBCO_TSBCO$contrast$Diff, na.rm = TRUE)),
  p_value = c(contrast_TSBCM_TSM$t_test$p.value,
              contrast_IBCM_TSBCM$t_test$p.value,
              contrast_TSBCO_TSO$t_test$p.value,
              contrast_IBCO_TSBCO$t_test$p.value)
)

print(ndvi_results)

# Boxplots for each NDVI contrast
plot_contrast_box <- function(contrast, trt1, trt2){
  contrast$contrast %>%
    tidyr::pivot_longer(cols = c(trt1, trt2), names_to = "Treatment", values_to = "Delta") %>%
    ggplot(aes(x = Treatment, y = Delta, fill = Treatment)) +
    geom_boxplot() +
    labs(title = paste("Δ NDVI:", trt1, "vs", trt2), y = "Δ NDVI") +
    theme(legend.position = "none")
}

plot_contrast_box(contrast_TSBCM_TSM, "TS_BCM", "TS_M")
plot_contrast_box(contrast_IBCM_TSBCM, "I_BCM", "TS_BCM")
plot_contrast_box(contrast_TSBCO_TSO, "TS_BCO", "TS_O")
plot_contrast_box(contrast_IBCO_TSBCO, "I_BCO", "TS_BCO")
# I think the full plot is considerably easier to see the full analysis... I would keep just that and maybe put the contrast ones in the appendix.
# -----------------------------
# 5. Per-Tree Analyses
# -----------------------------
# 5.1 Matured Pods
mod_mature <- glmmTMB::glmmTMB(Matured.pod ~ Treatment*Season + Diameter + (1|Farm) + (1|TreeUID),
                               data = df, family = poisson) # This works without error,
# but the individual TreeUID is still wrong! This invalidates your results!!!

summary(mod_mature)
emmeans_mature <- emmeans(mod_mature, pairwise ~ Treatment|Season)
plot(emmeans_mature$emmeans)

# Predictive Plot: Biochar vs No Biochar with Seasons
# -----------------------------

# Make sure Biochar grouping exists in your data
df <- df %>%
  mutate(
    Biochar = ifelse(grepl("BC", Treatment), "With Biochar", "Without Biochar")
  )

# Fit Negative Binomial Model with Season and Biochar

# This model would be correct, but MUST be run on the corrected data (i.e., corrected TreeUID)
mod_mature <- glmmTMB::glmmTMB(Matured.pod ~ Treatment*Season + Diameter + (1|Farm) + (1|TreeUID),
                               data = df,
                               family = nbinom2)
# Use Negative Binomial for counts

# 1. Recode Treatment into 2 groups
df$Treatment2 <- ifelse(grepl("BC", df$Treatment), "With Biochar", "Without Biochar")
df$Treatment2 <- factor(df$Treatment2, levels = c("Without Biochar", "With Biochar"))

# 2. Refit model with Treatment2
mod_mature <- glmmTMB::glmmTMB(Matured.pod ~ Treatment2*Season + Diameter + (1|Farm) + (1|TreeUID),
                               data = df,
                               family = nbinom2)

# 3. Build prediction grid
pred_grid <- expand.grid(
  Diameter = seq(min(df$Diameter, na.rm = TRUE),
                 max(df$Diameter, na.rm = TRUE),
                 length.out = 100),
  Treatment2 = levels(df$Treatment2),
  Season = levels(df$Season)
)

# 4. Get predictions
pred_grid$Predicted <- predict(mod_mature,
                               newdata = pred_grid,
                               type = "response",
                               re.form = NA)

# 5. Plot
ggplot() +
  geom_point(data = df, aes(x = Diameter, y = Matured.pod, color = Treatment2), alpha = 0.5) +
  geom_line(data = pred_grid, aes(x = Diameter, y = Predicted, color = Treatment2), size = 1.3) +
  facet_wrap(~Season) +
  labs(x = "Tree Diameter (cm)", y = "Predicted Matured Pod Count", color = "Treatment") +
  scale_color_manual(values = c("Without Biochar" = "red", "With Biochar" = "blue")) +
  theme_bw()


# ANOVA table
anova_mature <- car::Anova(mod_mature, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_mature)             # Convert to data frame
write.csv(anova_df, "Matured_pods_anova_results.csv", row.names = FALSE)
# ---> Make sure you can do a parametric comparison here! Same for all ANOVA and t-tests!

# EMMeans results
emmeans_mature <- emmeans(mod_mature, pairwise ~ Treatment2|Season)

# Save estimated marginal means
emmeans_df <- as.data.frame(emmeans_mature$emmeans)
write.csv(emmeans_df, "Matured_pods_emmeans_results.csv", row.names = FALSE)

# Save pairwise contrasts if needed
contrast_df <- as.data.frame(emmeans_mature$contrasts)
write.csv(contrast_df, "Matured_pods_pairwise_contrasts.csv", row.names = FALSE)



# 5.2 Young Pods
mod_young <- glmer.nb(No..of.young.pods ~ Treatment*Season + Diameter + (1|Farm) + (1|TreeUID), data = df)
summary(mod_young)
emmeans_young <- emmeans(mod_young, pairwise ~ Treatment|Season)
plot(emmeans_young$emmeans)

# Predictive Plot: Biochar vs No Biochar with Seasons
# -----------------------------

# 1. Recode Treatment into 2 groups
df$Treatment2 <- ifelse(grepl("BC", df$Treatment), "With Biochar", "Without Biochar")
df$Treatment2 <- factor(df$Treatment2, levels = c("Without Biochar", "With Biochar"))

# 2. Refit model with Treatment2
mod_young <- glmmTMB::glmmTMB(No..of.young.pods ~ Treatment2*Season + Diameter + (1|Farm) + (1|TreeUID),
                               data = df,
                               family = nbinom2)

# 3. Build prediction grid
pred_grid <- expand.grid(
  Diameter = seq(min(df$Diameter, na.rm = TRUE),
                 max(df$Diameter, na.rm = TRUE),
                 length.out = 100),
  Treatment2 = levels(df$Treatment2),
  Season = levels(df$Season)
)

# 4. Get predictions
pred_grid$Predicted <- predict(mod_young,
                               newdata = pred_grid,
                               type = "response",
                               re.form = NA)

# 5. Plot
ggplot() +
  geom_point(data = df, aes(x = Diameter, y = No..of.young.pods, color = Treatment2), alpha = 0.5) +
  geom_line(data = pred_grid, aes(x = Diameter, y = Predicted, color = Treatment2), size = 1.3) +
  facet_wrap(~Season) +
  labs(x = "Tree Diameter (cm)", y = "Predicted Young Pod Count", color = "Treatment") +
  scale_color_manual(values = c("Without Biochar" = "red", "With Biochar" = "blue")) +
  theme_bw()



# ANOVA table
anova_young <- car::Anova(mod_young, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_young)             # Convert to data frame
write.csv(anova_df, "Young_pods_anova_results.csv", row.names = FALSE)

# EMMeans results
emmeans_young <- emmeans(mod_young, pairwise ~ Treatment2|Season)

# Save estimated marginal means
emmeans_df <- as.data.frame(emmeans_young$emmeans)
write.csv(emmeans_df, "Young_pods_emmeans_results.csv", row.names = FALSE)

# Save pairwise contrasts if needed
contrast_df <- as.data.frame(emmeans_young$contrasts)
write.csv(contrast_df, "Young_pods_pairwise_contrasts.csv", row.names = FALSE)


# 5.3 Canopy Health
mod_canopy <- lmer(as.numeric(Canopy.health) ~ Treatment*Season + (1|Farm), data = df)
summary(mod_canopy)
# ---> Where is your TreeUID??? You can't ignore it. The canopy health from the
# wet season of Tree 1 is NOT INDEPENDENT of the canopy health of Tree 1 in the
# dry season!

# ANOVA table
anova_canopy <- car::Anova(mod_canopy, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_young)             # Convert to data frame
write.csv(anova_df, "Canopy_anova_results.csv", row.names = FALSE)


# 5.4 Flower Intensity
mod_flower <- lmer(as.numeric(Flower.Intensity) ~ Treatment*Season + (1|Farm), data = df)
summary(mod_flower)
# ---> Same here! Don't ignore TreeUID!

# ANOVA table
anova_flower <- car::Anova(mod_flower, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_flower)             # Convert to data frame
write.csv(anova_df, "Flower_anova_results.csv", row.names = FALSE)


# -----------------------------
# 6. Soil Δ Analysis
# -----------------------------
soil_df <- df %>%
  group_by(Farm, Treatment, Season) %>%
  summarise(across(starts_with("Mean."), \(x) mean(x, na.rm = TRUE)), .groups="drop") %>%
  tidyr::pivot_wider(names_from=Season, values_from=starts_with("Mean.")) %>%
  mutate(across(ends_with("Dry"), ~ . - get(sub("Dry","Wet",cur_column())), .names="Delta_{.col}"))

head(soil_df)


# 1. Compute Δ soil variables (Dry − Wet)
# -----------------------------
soil_delta <- df %>%
  group_by(Farm, Treatment, Season) %>%
  summarise(across(c(Mean.pH, Mean.N...., Mean.P..mg.Kg., Mean.K..mg.Kg.,
                     Mean..Ca.mg.Kg., Mean.Mg.mg.Kg., Mean.O.C...., Mean.O.M....),
                   mean, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Season,
                     values_from = c(Mean.pH, Mean.N...., Mean.P..mg.Kg., Mean.K..mg.Kg.,
                                     Mean..Ca.mg.Kg., Mean.Mg.mg.Kg., Mean.O.C...., Mean.O.M....)) %>%
  mutate(
    Delta_pH = `Mean.pH_Dry` - `Mean.pH_Wet`,
    Delta_N = `Mean.N...._Dry` - `Mean.N...._Wet`,
    Delta_P = `Mean.P..mg.Kg._Dry` - `Mean.P..mg.Kg._Wet`,
    Delta_K = `Mean.K..mg.Kg._Dry` - `Mean.K..mg.Kg._Wet`,
    Delta_Ca = `Mean..Ca.mg.Kg._Dry` - `Mean..Ca.mg.Kg._Wet`,
    Delta_Mg = `Mean.Mg.mg.Kg._Dry` - `Mean.Mg.mg.Kg._Wet`,
    Delta_OC = `Mean.O.C...._Dry` - `Mean.O.C...._Wet`,
    Delta_OM = `Mean.O.M...._Dry` - `Mean.O.M...._Wet`,
    Biochar = ifelse(grepl("BCM|BCO|BC", Treatment), "With_Biochar", "Without_Biochar")
  )

# -----------------------------
# 2. Select Δ variables
# -----------------------------
delta_vars <- dplyr::select(soil_delta, starts_with("Delta_"))

# ---> This is great so far! But... (see below)


# The formula below specifies that all the Delta variables are the response
# And Treatment and Farm are the predictors
manova_model <- manova(cbind(Delta_pH, Delta_N, Delta_P, Delta_K, Delta_Ca, Delta_Mg, Delta_OC, Delta_OM) ~ Treatment + Farm,
                       data = soil_delta)
summary(manova_model, test = "Pillai") # Use Pillai's Trace for robustness
# This would probably be a more adequate way to handle this.
summary.aov(manova_model)


# We create separate ANOVA models for each significant response variable.
# This makes the post-hoc analysis cleaner and more explicit.

# ANOVA for Delta_N
aov_N <- aov(Delta_N ~ Treatment + Farm, data = soil_delta)

# ANOVA for Delta_Mg
aov_Mg <- aov(Delta_Mg ~ Treatment + Farm, data = soil_delta)

# ANOVA for Delta_OC
aov_OC <- aov(Delta_OC ~ Treatment + Farm, data = soil_delta)

# --- Post-Hoc Test for Delta_N ---
tukey_N <- TukeyHSD(aov_N, "Treatment")
print(tukey_N)

# --- Post-Hoc Test for Delta_Mg ---
tukey_Mg <- TukeyHSD(aov_Mg, "Treatment")
print(tukey_Mg)

# --- Post-Hoc Test for Delta_OC ---
tukey_OC <- TukeyHSD(aov_OC, "Treatment")
print(tukey_OC)

# -----------------------------
# 7. Save Results
# -----------------------------
write.csv(ndvi_results, "output/ndvi_results.csv", row.names = FALSE)
write.csv(soil_df, "output/soil_delta.csv", row.names = FALSE)

# Optional: Save NDVI plots
ggsave("figures/ndvi_boxplot.png", width = 6, height = 4)

