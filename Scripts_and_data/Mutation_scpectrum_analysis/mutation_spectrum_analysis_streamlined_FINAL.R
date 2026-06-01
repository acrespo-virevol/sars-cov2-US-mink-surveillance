# =============================================================================
# MUTATION SPECTRUM ANALYSIS PIPELINE
# SARS-CoV-2 mink study — AY.39 lineage
# Single mink outbreak with known TMRCA
#
# INPUT FILES:
#   nextclade_mink.tsv  / nextclade_human.tsv  / nextclade_deer.tsv
#   meta_mink.csv       / meta_human.csv        / meta_deer.csv
#     (two columns each: seqName, samplingDate YYYY-MM-DD)
#   deer_metadata.csv
#     (columns: seqName, deer_outbreak_id, deer_origin)
#     deer_origin = "mink_derived" or "human_derived"
#
# OUTPUT FILES (written to OUTPUT_DIR):
#   CT_temporal_slopes_by_host.pdf
#   all_substitution_slopes_by_host.pdf
#   all_substitution_slopes_by_host.csv
#   CT_slope_by_host.csv
#   spectrum_boxplots.pdf
#   pca_mutation_spectrum.pdf
#   pca_deer_heterogeneity.pdf
#   mink_derived_deer_pca.pdf
#   mink_derived_deer_centroid_distances.csv
#   betadisper_diagnostics.pdf
#   betadisper_results.csv
#   pairwise_adonis_results.csv
#   dose_response_tmrca.pdf
#   dose_response_sensitivity.pdf
#   dose_response_data.csv
#   spectrum_proportions.csv
#   pca_scores.csv
#   mutation_spectrum_results.rds
# =============================================================================


# =============================================================================
# 0. SETUP
# =============================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  tidyverse,
  vegan,
  ggplot2,
  ggrepel,
  patchwork,
  lubridate,
  broom,
  RColorBrewer,
  scales
)

# ── USER CONFIGURATION ────────────────────────────────────────────────────────

MINK_FILE      <- "nextclade_mink.tsv"
HUMAN_FILE     <- "nextclade_human.tsv"
DEER_FILE      <- "nextclade_deer.tsv"

META_MINK      <- "meta_mink.csv"
META_HUMAN     <- "meta_human.csv"
META_DEER      <- "meta_deer.csv"

DEER_META_FILE <- "deer_metadata.csv"

OUTPUT_DIR        <- "mutation_spectrum_output_streamlined_FINAL"
MINK_TMRCA        <- as.Date("2022-05-29")   # <-- EDIT THIS
MIN_SUBSTITUTIONS <- 15
PANDEMIC_START    <- as.Date("2019-12-01")

HOST_COLORS <- c(
  "human"             = "#999999",
  "mink"              = "#4363d8",
  "deer"              = "#996600",
  "deer_mink_derived" = "#762A83"
)

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)
out_path <- function(f) file.path(OUTPUT_DIR, f)

cat("TMRCA:", format(MINK_TMRCA), "\n")


# =============================================================================
# 1. LOAD DATA
# =============================================================================

cat("\n=== STEP 1: Loading data ===\n")

load_nextclade <- function(filepath, host_label) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  df <- read_tsv(filepath, show_col_types = FALSE)
  
  name_col      <- intersect(c("seqName", "name", "sequence_name"),            names(df))[1]
  sub_col       <- intersect(c("substitutions", "nucSubstitutions"),            names(df))[1]
  total_sub_col <- intersect(c("totalSubstitutions", "totalNucSubstitutions"), names(df))[1]
  
  if (any(is.na(c(name_col, sub_col, total_sub_col)))) {
    cat("WARNING: expected columns not found in", filepath, "\n")
    cat("Available:", paste(names(df), collapse = ", "), "\n")
  }
  
  df %>%
    rename(seqName            = !!name_col,
           substitutions      = !!sub_col,
           totalSubstitutions = !!total_sub_col) %>%
    select(seqName, substitutions, totalSubstitutions) %>%
    mutate(host               = host_label,
           totalSubstitutions = as.integer(totalSubstitutions))
}

load_meta <- function(filepath) {
  if (!file.exists(filepath)) stop(paste("File not found:", filepath))
  read_csv(filepath, show_col_types = FALSE) %>%
    mutate(samplingDate     = as.Date(samplingDate),
           days_since_start = as.numeric(samplingDate - PANDEMIC_START)) %>%
    select(seqName, samplingDate, days_since_start)
}

mink_raw  <- load_nextclade(MINK_FILE,  "mink")  %>% left_join(load_meta(META_MINK),  by = "seqName")
human_raw <- load_nextclade(HUMAN_FILE, "human") %>% left_join(load_meta(META_HUMAN), by = "seqName")
deer_raw  <- load_nextclade(DEER_FILE,  "deer")  %>% left_join(load_meta(META_DEER),  by = "seqName")

# Remove sequences without dates
for (nm in c("mink_raw", "human_raw", "deer_raw")) {
  df  <- get(nm)
  bad <- sum(is.na(df$samplingDate))
  if (bad > 0) cat("WARNING:", bad, nm, "sequences missing dates — excluded\n")
  assign(nm, df %>% filter(!is.na(samplingDate)))
}

# ── Deer metadata ─────────────────────────────────────────────────────────────

deer_meta <- read_csv(DEER_META_FILE, show_col_types = FALSE) %>%
  mutate(deer_origin = coalesce(deer_origin, "human_derived"))

deer_raw <- deer_raw %>% left_join(deer_meta, by = "seqName")

# Separate mink-derived deer — excluded from all analyses, projected post-hoc
mink_derived_deer <- deer_raw %>%
  filter(deer_origin == "mink_derived") %>%
  mutate(host = "deer_mink_derived")

deer_raw_clean <- deer_raw %>%
  filter(deer_origin != "mink_derived")

cat("Mink-derived deer separated:", nrow(mink_derived_deer), "\n")

# ── Combine all groups ────────────────────────────────────────────────────────

all_seqs <- bind_rows(mink_raw, human_raw, deer_raw_clean, mink_derived_deer) %>%
  mutate(host = factor(host,
                       levels = c("human", "mink", "deer", "deer_mink_derived")))

cat("\nSequences loaded:\n")
count(all_seqs, host) %>% print()

cat("\nDate ranges:\n")
all_seqs %>%
  group_by(host) %>%
  summarise(earliest = min(samplingDate), latest = max(samplingDate),
            n = n(), .groups = "drop") %>%
  print()


# =============================================================================
# 2. PARSE SUBSTITUTIONS
# =============================================================================

cat("\n=== STEP 2: Parsing substitution vectors ===\n")

SUB_TYPES <- c("C>T","C>A","C>G","T>C","T>A","T>G",
               "G>A","G>T","G>C","A>G","A>T","A>C")

parse_substitutions <- function(sub_string) {
  counts <- setNames(rep(0L, 12), SUB_TYPES)
  if (is.na(sub_string) || nchar(trimws(sub_string)) == 0) return(counts)
  subs  <- trimws(str_split(sub_string, ",")[[1]])
  subs  <- subs[nchar(subs) > 0]
  if (length(subs) == 0) return(counts)
  calls <- paste0(str_sub(subs, 1, 1), ">", str_sub(subs, -1, -1))
  for (s in calls[calls %in% SUB_TYPES]) counts[s] <- counts[s] + 1L
  return(counts)
}

cat("Parsing", nrow(all_seqs), "sequences...\n")
sub_list <- vector("list", nrow(all_seqs))
for (i in seq_len(nrow(all_seqs))) {
  if (i %% 500 == 0) cat("  ", i, "of", nrow(all_seqs), "\n")
  sub_list[[i]] <- parse_substitutions(all_seqs$substitutions[i])
}
sub_matrix           <- do.call(rbind, sub_list)
rownames(sub_matrix) <- all_seqs$seqName
all_seqs             <- all_seqs %>% bind_cols(as_tibble(sub_matrix))
cat("Done.\n")


# =============================================================================
# 3. FILTER AND NORMALIZE TO PROPORTIONS
# =============================================================================

cat("\n=== STEP 3: Filter and normalize ===\n")

n_before <- nrow(all_seqs)
all_seqs <- all_seqs %>%
  filter(!is.na(totalSubstitutions), totalSubstitutions >= MIN_SUBSTITUTIONS)
cat("Removed", n_before - nrow(all_seqs),
    "sequences with <", MIN_SUBSTITUTIONS, "substitutions\n")
cat("Remaining:\n")
count(all_seqs, host) %>% print()

count_mat   <- as.matrix(all_seqs[, SUB_TYPES])
prop_mat    <- sweep(count_mat, 1, rowSums(count_mat), "/")

spectrum_df <- all_seqs %>%
  select(seqName, host, samplingDate, days_since_start,
         deer_outbreak_id, deer_origin) %>%
  bind_cols(as_tibble(prop_mat))

# Convenience subsets
main_df    <- spectrum_df %>% filter(host != "deer_mink_derived")
derived_df <- spectrum_df %>% filter(host == "deer_mink_derived")

cat("\nC>T proportion by host:\n")
spectrum_df %>%
  group_by(host) %>%
  summarise(n = n(), mean = round(mean(`C>T`), 4),
            sd = round(sd(`C>T`), 4), .groups = "drop") %>%
  print()


# =============================================================================
# 4. TEMPORAL SLOPE ANALYSIS
# =============================================================================
# No temporal correction is applied — PCA uses raw proportions.
# This section characterizes the host-specific temporal dynamics
# of C>T (and all 12 types) as a biological finding in its own right.

cat("\n=== STEP 4: Temporal slope analysis ===\n")

# ── C>T interaction model ─────────────────────────────────────────────────────

cat("C>T host × time interaction:\n")
interaction_model <- lm(
  `C>T` ~ days_since_start * host,
  data = main_df
)
summary(interaction_model) %>% print()

# Per-host slopes
slopes_by_host <- main_df %>%
  group_by(host) %>%
  group_map(~ {
    m <- lm(`C>T` ~ days_since_start, data = .x)
    tibble(host      = .y$host,
           slope     = coef(m)[2],
           slope_se  = summary(m)$coefficients[2, 2],
           p_value   = summary(m)$coefficients[2, 4],
           r_squared = summary(m)$r.squared)
  }) %>%
  bind_rows()

cat("\nC>T slopes by host:\n")
print(slopes_by_host)
write_csv(slopes_by_host, out_path("CT_slope_by_host.csv"))

# Percentage decline in mink
mink_ct_model <- lm(`C>T` ~ days_since_start,
                    data = spectrum_df %>% filter(host == "mink"))
ct_at_start   <- coef(mink_ct_model)[1]
mink_slope    <- slopes_by_host %>% filter(host == "mink") %>% pull(slope)
outbreak_days <- as.numeric(
  max(spectrum_df$samplingDate[spectrum_df$host == "mink"]) -
    min(spectrum_df$samplingDate[spectrum_df$host == "mink"])
)
pct_decline <- (mink_slope * outbreak_days / ct_at_start) * 100

cat("\nMink C>T percentage change over outbreak:\n")
cat("  Outbreak span:", outbreak_days, "days\n")
cat("  Starting C>T at outbreak onset:", round(ct_at_start, 4), "\n")
cat("  Total decline:", round(mink_slope * outbreak_days, 4), "\n")
cat("  Percentage change:", round(pct_decline, 1), "%\n")

# ── C>T slopes plot — main text Figure 6 ─────────────────────────────────────

# Combine main groups + mink-derived deer for plotting
ct_plot_df <- bind_rows(main_df, derived_df) %>%
  mutate(host = factor(host,
                       levels = c("human", "mink", "deer", "deer_mink_derived")))

# Per-group slope stats for all four groups (used for in-plot annotations)
slopes_ct_plot <- ct_plot_df %>%
  group_by(host) %>%
  group_map(~ {
    m <- lm(`C>T` ~ days_since_start, data = .x)
    tibble(
      host      = .y$host,
      r_squared = summary(m)$r.squared,
      p_value   = summary(m)$coefficients[2, 4],
      max_date  = max(.x$samplingDate),
      y_pred    = predict(m, newdata = data.frame(
        days_since_start = max(.x$days_since_start)))
    )
  }) %>%
  bind_rows() %>%
  mutate(
    host  = factor(host,
                   levels = c("human", "mink", "deer", "deer_mink_derived")),
    label = paste0("R² = ", round(r_squared, 3),
                   "\np = ", formatC(p_value, format = "e", digits = 2))
  )

p_ct_slopes <- ggplot(
  ct_plot_df,
  aes(x = samplingDate, y = `C>T`, color = host)
) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  geom_label(
    data        = slopes_ct_plot,
    aes(x       = max_date, y = y_pred,
        label   = label, color = host),
    size        = 3.2, hjust = 1.08, vjust = 0.5,
    fill        = "white", alpha = 0.85, label.size = 0.3,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = HOST_COLORS,
    labels = c("human"             = "Human",
               "mink"              = "Mink",
               "deer"              = "Deer",
               "deer_mink_derived" = "Deer (mink-derived)")
  ) +
  scale_x_date(date_labels = "%b %Y") +
  labs(x = "Collection date", y = "C>T proportion", color = "Host") +
  theme_bw(base_size = 16) +
  theme(legend.position = "right")

ggsave(out_path("CT_temporal_slopes_by_host.pdf"),
       p_ct_slopes, width = 10, height = 7)
cat("Saved: CT_temporal_slopes_by_host.pdf\n")

# ── All 12 substitution types — supplementary ─────────────────────────────────

all_slopes <- map_dfr(SUB_TYPES, function(st) {
  map_dfr(c("human", "mink", "deer"), function(h) {
    sub_data <- spectrum_df %>% filter(host == h)
    m        <- lm(as.formula(paste0("`", st, "` ~ days_since_start")),
                   data = sub_data)
    tibble(
      substitution_type = st,
      host              = h,
      slope             = coef(m)[2],
      p_value           = summary(m)$coefficients[2, 4],
      r_squared         = summary(m)$r.squared,
      sig_label         = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE            ~ ""
      )
    )
  })
})

write_csv(all_slopes, out_path("all_substitution_slopes_by_host.csv"))

cat("\nSignificant slopes (p < 0.05):\n")
all_slopes %>% filter(p_value < 0.05) %>% arrange(p_value) %>% print()

# Label positions for asterisks — placed at end of significant regression lines
label_positions <- all_slopes %>%
  filter(sig_label != "") %>%
  left_join(
    spectrum_df %>%
      filter(host != "deer_mink_derived") %>%
      group_by(host) %>%
      summarise(max_date = max(samplingDate),
                max_days = max(days_since_start),
                .groups  = "drop"),
    by = "host"
  ) %>%
  left_join(
    map_dfr(SUB_TYPES, function(st) {
      map_dfr(c("human", "mink", "deer"), function(h) {
        sub_data <- spectrum_df %>% filter(host == h)
        m  <- lm(as.formula(paste0("`", st, "` ~ days_since_start")),
                 data = sub_data)
        tibble(substitution_type = st, host = h,
               label_y = predict(m, newdata = data.frame(
                 days_since_start = max(sub_data$days_since_start))))
      })
    }),
    by = c("substitution_type", "host")
  ) %>%
  mutate(label_y_nudge = label_y * 1.10)

spectrum_long_time <- spectrum_df %>%
  filter(host != "deer_mink_derived") %>%
  select(seqName, host, samplingDate, all_of(SUB_TYPES)) %>%
  pivot_longer(all_of(SUB_TYPES),
               names_to  = "substitution_type",
               values_to = "proportion") %>%
  mutate(host = factor(host, levels = c("human", "mink", "deer")))

p_all_slopes <- ggplot(spectrum_long_time,
                       aes(x = samplingDate, y = proportion,
                           color = host)) +
  geom_point(alpha = 0.12, size = 1) +
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 0.85, alpha = 0.12) +
  geom_text(
    data        = label_positions,
    aes(x       = max_date, y = label_y_nudge,
        label   = sig_label, color = host),
    size        = 3.5, fontface = "bold",
    inherit.aes = FALSE, show.legend = FALSE
  ) +
  scale_color_manual(values = HOST_COLORS) +
  scale_x_date(date_labels = "%b %Y",
               guide = guide_axis(angle = 35)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ substitution_type, scales = "free_y", ncol = 4) +
  labs(#subtitle = "* p<0.05  ** p<0.01  *** p<0.001 at line end",
    x = "Collection date",
    y = "Proportion of substitutions",
    color = "Host") +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        strip.text       = element_text(face = "bold", size = 9))

ggsave(out_path("all_substitution_slopes_by_host.pdf"),
       p_all_slopes, width = 9, height = 10)
cat("Saved: all_substitution_slopes_by_host.pdf\n")


# =============================================================================
# 5. SPECTRUM BOXPLOTS
# =============================================================================

cat("\n=== STEP 5: Spectrum boxplots ===\n")

p_boxplots <- spectrum_df %>%
  filter(host %in% c("human", "mink", "deer")) %>%
  select(seqName, host, all_of(SUB_TYPES)) %>%
  pivot_longer(all_of(SUB_TYPES),
               names_to  = "substitution_type",
               values_to = "proportion") %>%
  ggplot(aes(x = host, y = proportion, fill = host)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, linewidth = 0.4) +
  scale_fill_manual(values = HOST_COLORS) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~ substitution_type, scales = "fixed", ncol = 4) +
  labs(x = "Host", y = "Proportion of substitutions", fill = "Host") +
  theme_bw(base_size = 10) +
  theme(legend.position  = "bottom",
        strip.background = element_rect(fill = "grey92"),
        axis.text.x      = element_text(angle = 30, hjust = 1))

ggsave(out_path("spectrum_boxplots.pdf"),
       p_boxplots, width = 8, height = 10)
cat("Saved: spectrum_boxplots.pdf\n")


# =============================================================================
# 6. ADONIS2 — MULTIVARIATE SPECTRUM COMPARISON
# =============================================================================
# Raw proportions with Bray-Curtis.
# Collection date included as covariate because mink and human sequences
# were collected across non-overlapping time periods.

cat("\n=== STEP 6: adonis2 ===\n")

prop_cols_main <- main_df %>% select(all_of(SUB_TYPES))
dist_matrix    <- vegdist(prop_cols_main, method = "bray")

set.seed(42)
adonis_overall <- adonis2(
  dist_matrix ~ days_since_start + host,
  data         = main_df,
  permutations = 9999,
  by           = "margin"
)
cat("=== Overall adonis2 ===\n")
print(adonis_overall)
host_row <- which(rownames(adonis_overall) == "host")
cat("Host R² =", round(adonis_overall$R2[host_row], 4),
    "(after accounting for collection date)\n")

run_pairwise <- function(g1, g2) {
  sub_meta  <- main_df %>% filter(host %in% c(g1, g2))
  sub_props <- spectrum_df %>%
    filter(host %in% c(g1, g2)) %>%
    select(all_of(SUB_TYPES))
  sub_dist <- vegdist(sub_props, method = "bray")
  set.seed(42)
  res <- adonis2(sub_dist ~ days_since_start + host,
                 data = sub_meta, permutations = 9999, by = "margin")
  hr  <- which(rownames(res) == "host")
  tibble(comparison = paste(g1, "vs", g2),
         n1      = sum(sub_meta$host == g1),
         n2      = sum(sub_meta$host == g2),
         R2      = round(res$R2[hr], 4),
         F_stat  = round(res$F[hr],  3),
         p_value = res$`Pr(>F)`[hr])
}

pairwise_results <- bind_rows(
  run_pairwise("mink", "human"),
  run_pairwise("deer", "human"),
  run_pairwise("mink", "deer")
) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "bonferroni"))

cat("\n=== Pairwise adonis2 (Bonferroni corrected) ===\n")
print(pairwise_results)
write_csv(pairwise_results, out_path("pairwise_adonis_results.csv"))
cat("Saved: pairwise_adonis_results.csv\n")


# =============================================================================
# 7. BETADISPER
# =============================================================================

cat("\n=== STEP 7: betadisper ===\n")

bd       <- betadisper(dist_matrix, group = main_df$host)
set.seed(42)
bd_perm  <- permutest(bd, permutations = 9999)
bd_tukey <- TukeyHSD(bd)

dispersion_summary <- tibble(
  host       = names(bd$group.distances),
  dispersion = as.numeric(bd$group.distances)
)
cat("Dispersions:\n")
print(dispersion_summary)
cat("Permutation test p =", round(bd_perm$tab$`Pr(>F)`[1], 4), "\n")
write_csv(dispersion_summary, out_path("betadisper_results.csv"))

bd_scores    <- as_tibble(bd$vectors) %>% mutate(host = main_df$host)
bd_centroids <- as_tibble(bd$centroids, rownames = "host")

p_bd <- ggplot() +
  geom_point(data = bd_scores,
             aes(x = PCoA1, y = PCoA2, color = host),
             alpha = 0.35, size = 1.4) +
  geom_point(data = bd_centroids,
             aes(x = PCoA1, y = PCoA2, color = host),
             size = 5, shape = 18) +
  geom_label_repel(data = bd_centroids,
                   aes(x = PCoA1, y = PCoA2, label = host, color = host),
                   size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = HOST_COLORS) +
  labs(subtitle = paste0("Permutation test p = ",
                         round(bd_perm$tab$`Pr(>F)`[1], 4)),
       x = "PCoA1", y = "PCoA2", color = "Host") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

p_bd_box <- ggplot(
  tibble(host = as.character(main_df$host), distance = bd$distances),
  aes(x = host, y = distance, fill = host)
) +
  geom_violin(alpha = 0.55, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.size = 0.5, color = "grey30") +
  scale_fill_manual(values = HOST_COLORS) +
  labs(x = "Host", y = "Distance to centroid") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave(out_path("betadisper_diagnostics.pdf"),
       p_bd / p_bd_box + plot_annotation(tag_levels = "A"),
       width = 10, height = 12)
cat("Saved: betadisper_diagnostics.pdf\n")


# =============================================================================
# 8. PCA — RAW PROPORTIONS
# =============================================================================

cat("\n=== STEP 8: PCA (raw proportions) ===\n")

pca_result    <- prcomp(prop_cols_main, scale. = TRUE, center = TRUE)
var_explained <- summary(pca_result)$importance[2, ] * 100
cumvar        <- cumsum(var_explained)
n_pcs_80      <- which(cumvar >= 80)[1]

cat("Variance explained:\n")
tibble(PC  = paste0("PC", 1:12),
       pct = round(var_explained[1:12], 2),
       cum = round(cumvar[1:12], 2)) %>% print()
cat("PCs for 80% variance:", n_pcs_80, "\n")

pca_scores_main <- as_tibble(pca_result$x) %>%
  bind_cols(main_df %>% select(seqName, host, samplingDate,
                               deer_outbreak_id))

# Project mink-derived deer post-hoc
if (nrow(derived_df) > 0) {
  mink_deer_projected <- predict(
    pca_result,
    newdata = derived_df %>% select(all_of(SUB_TYPES))
  ) %>%
    as_tibble() %>%
    bind_cols(derived_df %>% select(seqName, host, samplingDate,
                                    deer_outbreak_id))
  cat("Projected", nrow(mink_deer_projected),
      "mink-derived deer into PCA space\n")
} else {
  mink_deer_projected <- tibble()
}

group_centroids <- pca_scores_main %>%
  group_by(host) %>%
  summarise(centroid_PC1 = mean(PC1), centroid_PC2 = mean(PC2),
            .groups = "drop")

loading_scale <- min(max(abs(pca_scores_main$PC1)),
                     max(abs(pca_scores_main$PC2))) * 0.6

loadings_df <- as_tibble(pca_result$rotation, rownames = "substitution") %>%
  mutate(PC1_scaled = PC1 * loading_scale,
         PC2_scaled = PC2 * loading_scale,
         arrow_len  = sqrt(PC1^2 + PC2^2)) %>%
  arrange(desc(arrow_len)) %>%
  slice_head(n = 6)

# ── Main PCA biplot ────────────────────────────────────────────────────────────

p_pca <- ggplot() +
  geom_point(data  = pca_scores_main,
             aes(x = PC1, y = PC2, color = host),
             alpha = 0.4, size = 1.8) +
  stat_ellipse(data      = pca_scores_main,
               aes(x     = PC1, y = PC2, color = host),
               level = 0.95, linewidth = 0.9) +
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow     = arrow(length = unit(0.25, "cm"), type = "closed"),
               color     = "grey25", linewidth = 0.55) +
  geom_label_repel(data = loadings_df,
                   aes(x = PC1_scaled, y = PC2_scaled,
                       label = substitution),
                   size = 3.2, color = "grey15",
                   fill = alpha("white", 0.7),
                   label.size = 0.2, box.padding = 0.4) +
  scale_color_manual(values = HOST_COLORS) +
  labs(x     = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y     = paste0("PC2 (", round(var_explained[2], 1), "%)"),
       color = "Host") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")

p_pc1 <- ggplot(pca_scores_main, aes(x = host, y = PC1, fill = host)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.size = 0.5, color = "grey30") +
  scale_fill_manual(values = HOST_COLORS) +
  labs(x = "Host", y = "PC1") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_pc2 <- ggplot(pca_scores_main, aes(x = host, y = PC2, fill = host)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white",
               outlier.size = 0.5, color = "grey30") +
  scale_fill_manual(values = HOST_COLORS) +
  labs(x = "Host", y = "PC2") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

ggsave(out_path("pca_mutation_spectrum.pdf"),
       p_pca / (p_pc1 | p_pc2) +
         plot_layout(heights = c(2, 1)) +
         plot_annotation(tag_levels = "A"),
       width = 11, height = 13)
cat("Saved: pca_mutation_spectrum.pdf\n")

# ── Mink-derived deer overlay ──────────────────────────────────────────────────

if (nrow(mink_deer_projected) > 0) {
  
  centroid_distances <- mink_deer_projected %>%
    select(seqName, PC1, PC2) %>%
    cross_join(group_centroids) %>%
    mutate(dist = sqrt((PC1 - centroid_PC1)^2 +
                         (PC2 - centroid_PC2)^2)) %>%
    select(seqName, host, dist) %>%
    pivot_wider(names_from = host, values_from = dist,
                names_prefix = "dist_to_")
  
  cat("\nMink-derived deer distances to host centroids:\n")
  print(centroid_distances)
  write_csv(centroid_distances,
            out_path("mink_derived_deer_centroid_distances.csv"))
  
  p_mink_deer <- ggplot() +
    geom_point(data  = pca_scores_main,
               aes(x = PC1, y = PC2, color = host),
               alpha = 0.2, size = 1.2) +
    stat_ellipse(data      = pca_scores_main,
                 aes(x     = PC1, y = PC2, color = host),
                 level = 0.95, linewidth = 0.8, linetype = "dashed") +
    geom_point(data  = group_centroids,
               aes(x = centroid_PC1, y = centroid_PC2, color = host),
               size = 5, shape = 18) +
    geom_segment(
      data = mink_deer_projected %>%
        select(seqName, seq_PC1 = PC1, seq_PC2 = PC2) %>%
        cross_join(group_centroids),
      aes(x = seq_PC1, y = seq_PC2,
          xend = centroid_PC1, yend = centroid_PC2,
          color = host),
      linetype = "dotted", linewidth = 0.5, alpha = 0.6
    ) +
    geom_point(data  = mink_deer_projected,
               aes(x = PC1, y = PC2),
               color = HOST_COLORS["deer_mink_derived"],
               size = 4, shape = 17) +
    geom_label_repel(data = mink_deer_projected,
                     aes(x = PC1, y = PC2, label = seqName),
                     color = HOST_COLORS["deer_mink_derived"],
                     size = 3, box.padding = 0.5,
                     fill = alpha("white", 0.8)) +
    scale_color_manual(values = HOST_COLORS) +
    labs(x     = paste0("PC1 (", round(var_explained[1], 1), "%)"),
         y     = paste0("PC2 (", round(var_explained[2], 1), "%)"),
         color = "Host") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  
  ggsave(out_path("mink_derived_deer_pca.pdf"),
         p_mink_deer, width = 8, height = 6)
  cat("Saved: mink_derived_deer_pca.pdf\n")
}

# ── Deer outbreak heterogeneity ────────────────────────────────────────────────

deer_human_props <- spectrum_df %>%
  filter(host %in% c("human", "deer")) %>%
  select(all_of(SUB_TYPES))

deer_human_meta <- spectrum_df %>%
  filter(host %in% c("human", "deer"))

pca_deer_diag   <- prcomp(deer_human_props, scale. = TRUE, center = TRUE)
pca_deer_scores <- as_tibble(pca_deer_diag$x) %>%
  bind_cols(deer_human_meta %>% select(seqName, host, deer_outbreak_id))

n_outbreaks   <- n_distinct(na.omit(pca_deer_scores$deer_outbreak_id))
outbreak_ids  <- unique(na.omit(pca_deer_scores$deer_outbreak_id))
outbreak_cols <- setNames(
  colorRampPalette(brewer.pal(min(8, n_outbreaks), "Dark2"))(n_outbreaks),
  outbreak_ids
)

p_deer_het <- ggplot(pca_deer_scores, aes(x = PC1, y = PC2)) +
  geom_point(data  = pca_deer_scores %>% filter(host == "human"),
             color = "grey75", alpha = 0.2, size = 1.1) +
  geom_point(data  = pca_deer_scores %>% filter(host == "deer"),
             aes(color = deer_outbreak_id), size = 2.8, alpha = 0.85) +
  stat_ellipse(data = pca_deer_scores %>%
                 filter(host == "deer", !is.na(deer_outbreak_id)),
               aes(color = deer_outbreak_id),
               level = 0.80, linetype = "dashed") +
  scale_color_manual(values = outbreak_cols) +
  labs(x     = paste0("PC1 (",
                      round(summary(pca_deer_diag)$importance[2,1]*100, 1),
                      "%)"),
       y     = paste0("PC2 (",
                      round(summary(pca_deer_diag)$importance[2,2]*100, 1),
                      "%)"),
       color = "Deer outbreak") +
  theme_bw(base_size = 12) + theme(legend.position = "right")

ggsave(out_path("pca_deer_heterogeneity.pdf"),
       p_deer_het, width = 11, height = 8)
cat("Saved: pca_deer_heterogeneity.pdf\n")

write_csv(pca_scores_main, out_path("pca_scores.csv"))


# =============================================================================
# 9. DOSE-RESPONSE
# =============================================================================

cat("\n=== STEP 9: Dose-response ===\n")
cat("TMRCA:", format(MINK_TMRCA), "\n")
cat("Primary: PC1+PC2 | Sensitivity: PC1–PC", n_pcs_80,
    "(", round(cumvar[n_pcs_80], 1), "% variance)\n\n")

human_c_PC1 <- group_centroids %>% filter(host == "human") %>% pull(centroid_PC1)
human_c_PC2 <- group_centroids %>% filter(host == "human") %>% pull(centroid_PC2)

pc_cols <- paste0("PC", 1:n_pcs_80)
human_centroid_npcs <- pca_scores_main %>%
  filter(host == "human") %>%
  summarise(across(all_of(pc_cols), mean))

mink_dose <- pca_scores_main %>%
  filter(host == "mink") %>%
  mutate(
    days_since_tmrca = as.numeric(samplingDate - MINK_TMRCA),
    displacement_2pc = sqrt((PC1 - human_c_PC1)^2 +
                              (PC2 - human_c_PC2)^2),
    displacement_npc = sqrt(rowSums(
      across(all_of(pc_cols),
             ~ (. - human_centroid_npcs[[cur_column()]])^2)
    ))
  )

n_predating <- sum(mink_dose$days_since_tmrca < 0, na.rm = TRUE)
if (n_predating > 0) {
  cat("WARNING:", n_predating,
      "mink sequences predate TMRCA — excluded\n")
  mink_dose <- mink_dose %>% filter(days_since_tmrca >= 0)
}

cat("Mink sequences:", nrow(mink_dose), "\n")
cat("Circulation range:",
    min(mink_dose$days_since_tmrca), "–",
    max(mink_dose$days_since_tmrca), "days (",
    round(max(mink_dose$days_since_tmrca) / 30.4, 1), "months)\n")

# ── Project mink-derived deer into dose-response space ────────────────────────
# Uses the same human centroid and PCA space as mink; not included in regression.
deer_derived_dose <- mink_deer_projected %>%
  mutate(
    days_since_tmrca = as.numeric(samplingDate - MINK_TMRCA),
    displacement_2pc = sqrt((PC1 - human_c_PC1)^2 +
                              (PC2 - human_c_PC2)^2),
    displacement_npc = sqrt(rowSums(
      across(all_of(pc_cols),
             ~ (. - human_centroid_npcs[[cur_column()]])^2)
    ))
  ) %>%
  filter(days_since_tmrca >= 0)

cat("Mink-derived deer sequences overlaid on dose-response plot:",
    nrow(deer_derived_dose), "\n")

dose_2pc <- lm(displacement_2pc ~ days_since_tmrca, data = mink_dose)
dose_npc <- lm(displacement_npc ~ days_since_tmrca, data = mink_dose)

cat("\n=== PRIMARY (PC1+PC2) ===\n")
summary(dose_2pc) %>% print()

cat("\n=== SENSITIVITY (PC1–PC", n_pcs_80, ") ===\n")
cat("R² =", round(glance(dose_npc)$r.squared, 3),
    "  p =", formatC(tidy(dose_npc)$p.value[2], format = "e", digits = 2),
    "\n")

# Helper to build dose-response plot
make_dose_plot <- function(y_var, y_label, model, subtitle_extra = "",
                           overlay_df = NULL) {
  d_t <- tidy(model)
  d_g <- glance(model)
  
  ggplot(mink_dose, aes(x = days_since_tmrca, y = .data[[y_var]])) +
    geom_point(color = HOST_COLORS["mink"], alpha = 0.35, size = 2) +
    geom_smooth(method = "lm", color = "black",
                fill = "grey80", se = TRUE, linewidth = 1) +
    geom_vline(xintercept = 0, linetype = "dotted",
               color = "grey40", linewidth = 0.7) +
    annotate("text",
             x     = max(mink_dose$days_since_tmrca) * 0.02,
             y     = max(mink_dose[[y_var]]) * 0.97,
             label = "TMRCA", hjust = 0, size = 3.2, color = "grey40") +
    scale_x_continuous(
      name = "Days since TMRCA",
      sec.axis = sec_axis(~ . / 30.4,
                          name = "Months since TMRCA",
                          breaks = 0:24)
    ) +
    labs(
      subtitle = paste0(
        "R² = ", round(d_g$r.squared, 3),
        "  |  p = ", formatC(d_t$p.value[2], format = "e", digits = 2),
        "  |  n = ", nrow(mink_dose),
        subtitle_extra
      ),
      y = y_label
    ) +
    theme_bw(base_size = 16) -> p
  
  # Overlay mink-derived deer as a separate group (not in regression)
  if (!is.null(overlay_df) && nrow(overlay_df) > 0) {
    p <- p + geom_point(
      data        = overlay_df,
      aes(x       = days_since_tmrca, y = .data[[y_var]]),
      color       = HOST_COLORS["deer_mink_derived"],
      shape       = 17,   # triangle to distinguish from mink circles
      alpha       = 0.65, size = 2.5,
      inherit.aes = FALSE
    )
  }
  p
}

p_dose <- make_dose_plot(
  "displacement_2pc",
  "Euclidean distance from human centroid (PC1, PC2)",
  dose_2pc,
  "",
  overlay_df = deer_derived_dose
)

p_dose_sens <- make_dose_plot(
  "displacement_npc",
  paste0("Euclidean distance from human centroid (PC1–PC", n_pcs_80, ")"),
  dose_npc,
  paste0("  |  PC1–PC", n_pcs_80,
         " (", round(cumvar[n_pcs_80], 1), "% variance)"),
  overlay_df = deer_derived_dose
)

ggsave(out_path("dose_response_tmrca.pdf"),
       p_dose, width = 10, height = 7)
cat("Saved: dose_response_tmrca.pdf  (primary)\n")

ggsave(out_path("dose_response_sensitivity.pdf"),
       p_dose_sens, width = 10, height = 7)
cat("Saved: dose_response_sensitivity.pdf  (sensitivity)\n")

write_csv(
  mink_dose %>% select(seqName, samplingDate, days_since_tmrca,
                       PC1, PC2, all_of(pc_cols),
                       displacement_2pc, displacement_npc),
  out_path("dose_response_data.csv")
)
cat("Saved: dose_response_data.csv\n")


# =============================================================================
# 10. SAVE ALL OUTPUTS
# =============================================================================

cat("\n=== STEP 10: Saving results ===\n")

write_csv(spectrum_df, out_path("spectrum_proportions.csv"))

results_summary <- list(
  sample_sizes    = count(all_seqs, host),
  mink_tmrca      = MINK_TMRCA,
  slopes_CT       = slopes_by_host,
  pct_CT_decline  = pct_decline,
  adonis_overall  = adonis_overall,
  adonis_pairwise = pairwise_results,
  betadisper      = list(dispersions = dispersion_summary,
                         permtest    = bd_perm,
                         tukey       = bd_tukey),
  pca_variance    = tibble(PC  = paste0("PC", 1:12),
                           pct = round(var_explained[1:12], 3),
                           cum = round(cumvar[1:12], 3)),
  dose_response   = list(
    primary     = list(model = tidy(dose_2pc), fit = glance(dose_2pc)),
    sensitivity = list(model = tidy(dose_npc), fit = glance(dose_npc),
                       n_pcs = n_pcs_80,
                       var_pct = round(cumvar[n_pcs_80], 1))
  )
)

saveRDS(results_summary, out_path("mutation_spectrum_results.rds"))

cat("\nAll outputs saved to:", OUTPUT_DIR, "\n")
cat("\n=== PIPELINE COMPLETE ===\n")
cat("Output files:\n")
list.files(OUTPUT_DIR) %>% paste0("  ", .) %>% cat(sep = "\n")