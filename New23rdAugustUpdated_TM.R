# =========================================
# Biochar Cocoa Farm Analysis Script
# =========================================
# setwd("/Users/otissenalor/Documents/Otis Masters Thesis") # If you have a project, you don't need to set absolute paths.
# THIS AVOIDS PROBLEMS IN THE SECOND LINE OF CODE ALREADY. THAT IS THE BEAUTY OF THE Rproj. ;)
# -----------------------------
# 1. Setup
# -----------------------------
if(!require("Require")) install.packages("Require")
Require::Require(c("MASS","ggplot2","dplyr","lme4","car","emmeans","multcomp","tidyr"))

# Consistent ggplot theme
theme_set(theme_minimal(base_size = 14))

# -----------------------------
# 2. Load Data
# -----------------------------
df <- read.csv("csv2.csv", stringsAsFactors = TRUE)

# Create unique TreeID
df <- df %>%
  mutate(TreeUID = interaction(Farm, Tree.id, drop = TRUE))

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
contrast_TSBCO_TSO <- run_contrast(ndvi_df, "TS_BCO", "TS_O") # Here we have a problem
# You don't have TS_O for Farm 1. This invalidates the analysis. You should exclude FARM 1 from here.
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
mod_mature <- glmer.nb(Matured.pod ~ Treatment*Season + Diameter + (1|Farm) + (1|TreeUID), data = df)
# GOT AN ERROR HERE. DOUBLE CHECK! I think the reason is the fact that you assume a TreeUID for repeated measures
# (i.e., pre and post treatment) but in the end you only have one id per tree!
# > # 5.1 Matured Pods
#   > mod_mature <- glmer.nb(Matured.pod ~ Treatment*Season + Diameter + (1|Farm) + (1|TreeUID), data = df)
# Error in eval(mc, parent.frame(1L)) :
#   (maxstephalfit) PIRLS step-halvings failed to reduce deviance in pwrssUpdate
# But this works...
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
# --> Can you explain how you came to this model? What was it that you were trying to see?
# This is not the correct way of dealing with this.
model_biochar <- glm.nb(Matured.pod ~ Diameter * Biochar + Season, data = df)

# This is why this model was wrong:
# A glm is a "fixed effects only" model. It cannot handle random effects.
# This model completely ignores the most important features of the study design: It does not
# account for the blocking structure (Farm) or for the repeated measurements on the same tree (TreeUID).

# This model would be correct, but MUST be run on the corrected data (i.e., corrected TreeUID)
# mod_mature <- glmmTMB::glmmTMB(Matured.pod ~ Treatment*Season + Baseline_DBH + (1|Farm) + (1|TreeUID),
#                                data = your_new_corrected_df,
#                                family = nbinom2) # Use Negative Binomial for counts

# Create prediction grid
pred_grid <- expand.grid(
  Diameter = seq(min(df$Diameter, na.rm = TRUE), max(df$Diameter, na.rm = TRUE), length.out = 100),
  Biochar = c("With Biochar", "Without Biochar"),
  Season = unique(df$Season)
)

# Predict values
pred_grid$Predicted <- predict(model_biochar, newdata = pred_grid, type = "response")

# Plot faceted by Season
ggplot() +
  geom_point(data = df, aes(x = Diameter, y = Matured.pod, color = Biochar), alpha = 0.5) +
  geom_line(data = pred_grid, aes(x = Diameter, y = Predicted, color = Biochar), size = 1.3) +
  geom_vline(xintercept = 12, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 10, linetype = "dotted", color = "black") +
  facet_wrap(~Season) +
  labs(
    title = "Predicted Matured Pod Count by Diameter, Biochar, and Season",
    x = "Tree Diameter (cm)",
    y = "Predicted Matured Pod Count",
    color = "Biochar Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )


# ANOVA table
anova_mature <- car::Anova(mod_mature, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_mature)             # Convert to data frame
write.csv(anova_df, "Matured_pods_anova_results.csv", row.names = FALSE)
# ---> Make sure you can do a parametric comparison here! Same for all ANOVA and t-tests!

# EMMeans results
emmeans_mature <- emmeans(mod_mature, pairwise ~ Treatment|Season)

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

# Make sure Biochar grouping exists in your data
df <- df %>%
  mutate(
    Biochar = ifelse(grepl("BC", Treatment), "With Biochar", "Without Biochar")
  )

# Fit Negative Binomial Model with Season and Biochar
model_biochar <- glm.nb(No..of.young.pods ~ Diameter * Biochar + Season, data = df)
# ---> I can see what you are trying to do—create a simplified plot showing the 'Biochar vs. No Biochar' effect.
# This is a good idea. However, the way it's been done is statistically incorrect and we must fix it.
# But there are two issues here:
# 1. You are using the wrong type of model for the plot. The glm.nb model is a fixed-effects model.
# It completely ignores our experimental design—it doesn't account for the different farms or the
# fact that we measured the same trees twice. This leads to incorrect predictions, and the lines on
# your plot will not be accurate.
# 2. Your plot does not match your statistical test. Your ANOVA table is correctly calculated
# from your mixed-effects model, mod_young. But your plot is drawn from a different,
# incorrect model, model_biochar. In scientific reporting, your figures must always visualize
# the results from the exact same model you used for your hypothesis tests.
# You need to delete the model_biochar line entirely and generate the predictions
# for the plot directly from your main, correct model (mod_young).
# This ensures that your plot is an honest and accurate representation of your statistical findings.

#### IMPORTANT!!!! ####
# Please revise all the code where you are doing this and fix it!
#######################

# Create prediction grid
pred_grid <- expand.grid(
  Diameter = seq(min(df$Diameter, na.rm = TRUE), max(df$Diameter, na.rm = TRUE), length.out = 100),
  Biochar = c("With Biochar", "Without Biochar"),
  Season = unique(df$Season)
)

# Predict values
pred_grid$Predicted <- predict(model_biochar, newdata = pred_grid, type = "response")

# Plot faceted by Season
ggplot() +
  geom_point(data = df, aes(x = Diameter, y = No..of.young.pods, color = Biochar), alpha = 0.5) +
  geom_line(data = pred_grid, aes(x = Diameter, y = Predicted, color = Biochar), size = 1.3) +
  geom_vline(xintercept = 12, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 10, linetype = "dotted", color = "black") +
  facet_wrap(~Season) +
  labs(
    title = "Predicted Young Pod Count by Diameter, Biochar, and Season",
    x = "Tree Diameter (cm)",
    y = "Predicted Young Pod Count",
    color = "Biochar Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold")
  )

# ANOVA table
anova_young <- car::Anova(mod_young, type = 3)   # Type III ANOVA
anova_df <- broom::tidy(anova_young)             # Convert to data frame
write.csv(anova_df, "Young_pods_anova_results.csv", row.names = FALSE)

# EMMeans results
emmeans_young <- emmeans(mod_young, pairwise ~ Treatment|Season)

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

# -----------------------------
# 3. Fit LDA model
# -----------------------------
lda_model <- lda(Biochar ~ ., data = cbind(delta_vars, Biochar = soil_delta$Biochar))
# ---> This is wrong. The standard lda function in R has no mechanism to account for
# a blocking factor like Farm. It treats the 30 plots as 30 fully independent samples.
# This is statistically incorrect and violates the assumption of independence.
# This is a form of pseudoreplication. The model thinks it has more independent
# information than it really does because it's ignoring the grouping structure
# (Farm). This can make the separation between the groups look stronger or more
# "significant" than it truly is.

# The formula below specifies that all the Delta variables are the response
# And Treatment and Farm are the predictors
manova_model <- manova(cbind(Delta_pH, Delta_N, Delta_P, Delta_K, Delta_Ca, Delta_Mg, Delta_OC, Delta_OM) ~ Treatment + Farm,
                       data = soil_delta)
summary(manova_model, test = "Pillai") # Use Pillai's Trace for robustness
# This would probably be a more adequate way to handle this.

# -----------------------------
# 4. Cross-validated accuracy
# -----------------------------
lda_cv <- lda(Biochar ~ ., data = cbind(delta_vars, Biochar = soil_delta$Biochar), CV = TRUE)
cv_table <- table(Predicted = lda_cv$class, Actual = soil_delta$Biochar)
cv_accuracy <- sum(diag(cv_table)) / sum(cv_table)
cat("Cross-validated accuracy:", round(cv_accuracy*100, 1), "%\n")
print(cv_table)

# -----------------------------
# 5. Variable contributions (loadings)
# -----------------------------
loadings <- data.frame(lda_model$scaling)
loadings$Variable <- rownames(loadings)
loadings <- loadings %>%
  mutate(abs_LD1 = abs(LD1)) %>%
  arrange(desc(abs_LD1))

print(loadings)

# -----------------------------
# 6. LDA scores for plotting
# -----------------------------
lda_pred <- predict(lda_model)
lda_df <- data.frame(LD1 = lda_pred$x[,1],
                     Treatment = soil_delta$Treatment,
                     Biochar = soil_delta$Biochar)

# -----------------------------
# 7. Plot 1D LDA biplot
# -----------------------------
ggplot(lda_df, aes(x = LD1, y = 0, color = Biochar, shape = Treatment)) +
  geom_point(size = 3, alpha = 0.8, position = position_jitter(height = 0.05)) +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 15, 16, 17)) +  # assign shapes for 10 treatments
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Δ-LDA: Soil Variables (Dry − Wet) by Biochar",
       subtitle = paste("Cross-validated accuracy:", round(cv_accuracy*100,1), "%"),
       x = "LD1", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "right")

# -----------------------------
# 7. Save Results
# -----------------------------
write.csv(ndvi_results, "output/ndvi_results.csv", row.names = FALSE)
write.csv(soil_df, "output/soil_delta.csv", row.names = FALSE)

# Optional: Save NDVI plots
ggsave("figures/ndvi_boxplot.png", width = 6, height = 4)

