# ====================================================
# Primary Author: Kexin Gong
# Secondary Contributor: Eman Tahir
# Course: Software Tools (BINF 6210)

# Project: Comparing the BIN Composition of Subfamily Sciurinae between North America and Eurasia 

# Course: Software Tools (BINF 6210)
# ====================================================

## BACKGROUND:
# The subfamily Sciurinae (tree squirrels and flying squirrels) has a wide geographic range across North America and Eurasia. DNA barcode records from BOLD can be grouped into BINs (Barcode Index Numbers), which provide a standardized way to compare genetic lineages across regions.

## WHY THIS IS INTERESTING:
# Although Sciurinae species occur on both continents, their actual BIN diversity and overlap remain unclear. Differences in species ranges,  historical biogeography, and sampling effort may all shape how similar or distinct their BIN assemblages appear.

## RESEARCH QUESTION & HYPOTHESIS:
# How similar is the BIN composition of subfamily Sciurinae between North America and Eurasia?
# Hypothesis: Because of long-term geographic isolation, Sciurinae populations in North America and Eurasia will show distinct BIN assemblages with limited overlap.

## REPRODUCIBILITY NOTES:
# The script assumes the project root is open in an RStudio .Rproj so that 
# relative paths (data/, figs/) resolve correctly. 

## ======= 0: LOAD PACKAGES =======
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(readr)
library(stringr)

## ======= 1: READ DATA + QUICK CHECKS =======
# Read the BOLD result file and standardize the country column name.
dat <- read.delim("data/result.tsv", sep = "\t", header = TRUE, quote = "", check.names = FALSE)

# Replace the inconvenient column name used in the following code
dat <- dat %>% rename(country = `country/ocean`)

# Quick data-quality snapshot
cat("\n=== DATA QUALITY CHECK ===\n")
dat_summary <- dat %>%
  summarise(
    records         = n(),
    missing_country = sum(is.na(country)),
    missing_bin     = sum(is.na(bin_uri)),
    n_countries     = n_distinct(country),
    n_bins          = n_distinct(bin_uri)
  )
print(dat_summary)

## ======= 2: MAP COUNTRIES TO CONTINENTS =======
# Map countries from the BOLD export into two focal continents: North America and Eurasia.

# The following used country names based on the data file
north_america <- c(
  "Canada", "United States", "Mexico", "Panama", "Costa Rica"
)
eurasia <- c(
  "Austria", "China", "Denmark", "Georgia", "Germany", "India", "Italy",
  "Laos", "Lebanon", "Mongolia", "Norway", "Russia", "Slovakia", "South Korea",
  "Switzerland", "Taiwan", "Turkiye", "Vietnam"
)

# Assign each record to a continent
dat <- dat %>%
  mutate(continent = case_when(
    country %in% north_america ~ "North America",
    country %in% eurasia       ~ "Eurasia",
    TRUE ~ NA_character_
  ))

# Keep only the two continents of interest
dat2 <- dat %>% filter(continent %in% c("North America", "Eurasia"))

## ======= 3: BUILD BIN–BY–CONTINENT MATRICES =======
# Count how many records belong to each BIN in each continent, then:
#   (a) build a wide count matrix (BIN x continent)
#   (b) build a presence/absence matrix (1 = present, 0 = absent)

bin_counts <- dat2 %>%
  filter(!is.na(bin_uri)) %>%
  group_by(continent, bin_uri) %>%
  summarise(n = n(), .groups = "drop")

# Convert to wide format: one row per BIN, one column per continent
bin_wide <- bin_counts %>%
  tidyr::pivot_wider(names_from = continent, values_from = n, values_fill = 0)

# Presence/absence matrix (1 = present, 0 = absent)
# Scalable version using across(): works even if more continents are added later.
pa <- bin_wide %>%
  mutate(across(-bin_uri, ~ as.integer(.x > 0)))

## ======= 4: COMPUTE SIMILARITY METRICS =======
# Compare the two continents using:
# - Jaccard distance on presence/absence (composition only)
# - Bray-Curtis distance on counts (abundance-sensitive)
# Also count the number of shared BINs.

# Convert to 2×N matrix (rows = continents, columns = BINs)
mat <- as.matrix(t(pa[, -1]))
rownames(mat) <- colnames(pa[, -1])

# Compute binary Jaccard distance
dist_jac <- vegdist(mat, method = "jaccard", binary = TRUE)
print(dist_jac)

# Compute Bray–Curtis distance
mat_bray <- as.matrix(t(bin_wide %>% select(`North America`,`Eurasia`)))
dist_bray <- vegdist(mat_bray, method = "bray")

# Count shared BINs
shared_bins <- pa %>% filter(`North America` == 1 & `Eurasia` == 1) %>% pull(bin_uri)

message(sprintf("Unique BINs NA=%d, EUAS=%d, Shared=%d",
                sum(pa$`North America`), sum(pa$Eurasia), length(shared_bins)))

### --- Major Edit 1: Rarefied Bin Richness ---
# Motivation:
# Eurasia has more public records than North America, so simple BIN counts may be inflated just because of higher sampling effort. This individual-based rarefaction will allow us to compare expected BIN richness at a common sample size (the smaller total record count) without sampling bias.

# Total record counts per continent (from the count matrix used for Bray–Curtis)
total_records <- rowSums(mat_bray)
total_records

# Rarefy both continents down to the smaller sample size
rare_sample_size <- min(total_records)

# vegan::rarefy expects a community matrix (rows = samples, cols = species)
# and returns expected richness at 'sample' individuals.
rare_vals <- vegan::rarefy(mat_bray, sample = rare_sample_size)

# Build a tidy summary table of raw vs rarefied richness
bin_rich_raw <- pa %>%
  summarise(
    `North America` = sum(`North America`),
    `Eurasia`       = sum(Eurasia)
  ) %>%
  pivot_longer(everything(),
               names_to  = "continent",
               values_to = "raw_bin_richness")

rarefied_df <- bin_rich_raw %>%
  mutate(
    total_records      = as.integer(total_records[continent]),
    rarefaction_sample = rare_sample_size,
    rarefied_richness  = as.numeric(rare_vals[continent])
  )

cat("\n=== RAREFIED BIN RICHNESS (MAJOR EDIT 1) ===\n")
print(rarefied_df)

### --- Major Edit 2: Permutation test for shared BIN overlap ---
# Motivation:
# Jaccard and Bray-Curtis describe differences, but they do not tell us if the number of shared BINs is unusually high or low. Here, we shuffle BINs between continents (while keeping each continent’s total BIN count the same) to see what amount of overlap would happen by chance. This gives us a fair comparison for the observed overlap.

set.seed(6210)  # reproducible permutation results

# Observed overlap and BIN counts per continent
obs_shared <- length(shared_bins)
na_bins    <- sum(pa$`North America` == 1)
eu_bins    <- sum(pa$Eurasia == 1)
all_bins   <- pa$bin_uri

# Container for null draws (used later for Figure 5)
perm_df <- NULL

if (length(all_bins) > 0 && na_bins > 0 && eu_bins > 0) {
  B <- 5000L  # number of permutations
  
  # Null model: randomly assign BIN membership, keeping richness per continent
  sim_shared <- replicate(B, {
    rand_na <- sample(all_bins, na_bins, replace = FALSE)
    rand_eu <- sample(all_bins, eu_bins, replace = FALSE)
    length(intersect(rand_na, rand_eu))
  })
  
  # One-sided p-value: Pr(null <= observed overlap)  # lower-than-expected overlap
  p_val_low <- mean(sim_shared <= obs_shared)
  
  perm_summary <- tibble::tibble(
    observed_shared = obs_shared,
    null_mean       = mean(sim_shared),
    null_sd         = stats::sd(sim_shared),
    p_value_low     = p_val_low
  )

  
  cat("\n=== PERMUTATION TEST: SHARED BIN OVERLAP (MAJOR EDIT 2) ===\n")
  print(perm_summary)
  
  # Store null draws for visualization in the Figures section
  perm_df <- tibble::tibble(sim_shared = sim_shared)
} else {
  cat("\n=== PERMUTATION TEST SKIPPED: not enough BINs to simulate ===\n")
}

### --- Major Edit 3: BIN-level enrichment analysis (Fisher’s Exact Test + FDR) ---
# Motivation:
# Even if overall similarity is low, we want to identify WHICH specific BINs are disproportionately represented in one continent vs the other. Fisher’s Exact Test is ideal for categorical count data with small sample sizes.

# Total counts per continent (all BINs combined)
tot_na <- sum(bin_wide$`North America`)
tot_eu <- sum(bin_wide$Eurasia)

# Build BIN-level contingency tables and run Fisher’s test per BIN
bin_tests <- bin_wide %>%
  rowwise() %>%
  mutate(
    fisher_p = {
      mat <- matrix(
        c(`North America`,               # this BIN in NA
          Eurasia,                       # this BIN in Eurasia
          tot_na - `North America`,      # other BINs in NA
          tot_eu - Eurasia),             # other BINs in Eurasia
        nrow = 2,
        byrow = TRUE
      )
      fisher.test(mat)$p.value
    }
  ) %>%
  ungroup()

# Add FDR correction
bin_tests <- bin_tests %>%
  mutate(fdr = p.adjust(fisher_p, method = "BH"))

# Determine enrichment direction
bin_tests <- bin_tests %>%
  mutate(
    direction = case_when(
      fdr >= 0.05 ~ "NS",
      `North America` > Eurasia ~ "North America enriched",
      Eurasia > `North America` ~ "Eurasia enriched",
      TRUE ~ "NS"
    )
  )

cat("\n=== BIN-LEVEL ENRICHMENT RESULTS (MAJOR EDIT 3) ===\n")
print(bin_tests %>% select(bin_uri, `North America`, Eurasia, fisher_p, fdr, direction))

## ======= 5: FIGURES =======
# Figure 1: BIN richness per continent (bar plot)
bin_rich <- pa |>
  summarise(`North America` = sum(`North America`),
            `Eurasia` = sum(Eurasia)) |>
  pivot_longer(everything(), names_to = "continent", values_to = "unique_bins")

p1 <- ggplot(bin_rich, aes(continent, unique_bins, fill = continent)) +
  geom_col(width = 0.7) +
  labs(title = "Figure 1: Sciurinae BIN richness by continent",
       x = NULL, y = "Number of unique BINs") +
  scale_fill_manual(values = c("Eurasia" = "#FFC3C3", "North America" = "#FFE3A9"))+
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  ylim(0,15)

ggsave("figs/Fig1_BIN_richness_by_continent.png", p1, width = 6, height = 4, dpi = 300)

# Figure 2: Presence/absence heatmap (BIN × continent)
pa_long <- pa %>%
  pivot_longer(cols = c(`North America`, `Eurasia`),
               names_to = "continent", values_to = "present")

# Sort BINs by total occurrences to make the heatmap cleaner
keep_ids <- pa %>% mutate(tot = `North America` + Eurasia) %>%
  arrange(desc(tot)) %>% pull(bin_uri)
pa_long$bin_uri <- factor(pa_long$bin_uri, levels = keep_ids)

p2 <- ggplot(pa_long, aes(continent, bin_uri, fill = factor(present))) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("0" = "#eeeeee", "1" = "#B0CE88"), name = "Present") +
  labs(title = "Figure 2: Presence/absence of Sciurinae BINs across continents",
       x = NULL, y = "BIN") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())

ggsave("figs/Fig2_BIN_presence_heatmap.png", p2, width = 6, height = 6, dpi = 300)

# Figure 3: Number of public records by continent
rec_counts <- dat2 %>% count(continent)
p3 <- ggplot(rec_counts, aes(continent, n, fill = continent)) +
  geom_col(width = 0.7) +
  labs(title = "Figure 3: Number of public records by continent",
       x = NULL, y = "Records") +
  scale_fill_manual(values = c("Eurasia" = "#FFC3C3", "North America" = "#FFE3A9"))+
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank())

ggsave("figs/Fig3_public_records_by_continent.png", p3, width = 6, height = 4, dpi = 300)

# Figure 4: rarefied BIN richness (Optional visualization)
p4 <- ggplot(rarefied_df,
             aes(x = continent, y = rarefied_richness, fill = continent)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", total_records)),
            vjust = -0.4, size = 3.2) +
  labs(
    title = "Figure 4: Rarefied Sciurinae BIN richness by continent",
    subtitle = paste0("Expected richness at ", rare_sample_size,
                      " records per continent (individual-based rarefaction)"),
    x = NULL,
    y = "Rarefied BIN richness"
  ) +
  scale_fill_manual(values = c("Eurasia" = "#FFC3C3", "North America" = "#FFE3A9")) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid      = element_blank()
  )

ggsave("figs/Fig4_Rarefied_BIN_richness.png", p4, width = 6, height = 4, dpi = 300)

# Figure 5: Null distribution of shared BINs under random BIN assignment
if (!is.null(perm_df)) {
  p5 <- ggplot(perm_df, aes(x = sim_shared)) +
    geom_histogram(binwidth = 1, boundary = -0.5, colour = "white", fill = "#B0CE88") +
    geom_vline(xintercept = obs_shared, linewidth = 1.1) +
    labs(
      title = "Figure 5: Null distribution of shared BIN overlaps",
      subtitle = paste0("Observed overlap = ", obs_shared,
                        " (vertical line). Null model preserves BIN counts per continent."),
      x = "Shared BINs under null model",
      y = "Frequency"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank()
    )
  ggsave("figs/Fig5_Shared_BINs_null_distribution.png", p5, width = 6, height = 4, dpi = 300)
}

## ======= 6: SUMMARY OUTPUT =======
cat("\n=== SUMMARY ===\n")
cat("File:", "data/result.tsv", "\n")
cat("Records used:", nrow(dat2), "\n")
cat("Continents:\n"); print(table(dat2$continent))
cat("Unique BINs per continent:\n"); print(bin_rich)
cat("Jaccard distance (0 = identical, 1 = completely different):\n"); print(dist_jac)
cat("Bray–Curtis distance:\n"); print(as.matrix(dist_bray))
cat("Shared BIN IDs (first 10):\n"); print(head(shared_bins, 10))