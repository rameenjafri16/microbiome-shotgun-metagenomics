# 04_diversity.R
# Rameen Jafri
#
#
# Diversity analysis of shotgun metagenomics data
# Includes: taxonomic abundance visualization, alpha diversity, beta diversity
# Tools: phyloseq, vegan, ggplot2
#
# Run from: ~/BINF_DESKTOP/BINF6110/assignment03/
# Input: results/bracken/combined_abundance.txt

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


# 1. LOAD AND FORMAT DATA -----------------

# Load combined Bracken output
#setwd("~/BINF_DESKTOP/BINF6110/assignment03")
bracken <- read.table("results/bracken/combined_abundance.txt",
                      header = TRUE, sep = "\t", row.names = 1)

# Keep only fraction columns (remove raw count columns)
# Bracken combine output has columns: name, taxid, SRR_num, SRR_frac ...
# Keep only the *_frac columns for relative abundance
frac_cols <- grep("_frac", colnames(bracken))
counts_cols <- grep("_num", colnames(bracken))

# Use raw counts for DESeq2 (saved separately) and phyloseq
otu_counts <- bracken[, counts_cols]
colnames(otu_counts) <- c("SRR8146968", "SRR8146978", "SRR8146982",
                          "SRR8146938", "SRR8146935", "SRR8146936")

# Sample metadata
metadata <- data.frame(
  sample_id  = c("SRR8146968", "SRR8146978", "SRR8146982",
                 "SRR8146938", "SRR8146935", "SRR8146936"),
  diet_group = c("Vegan", "Vegan", "Vegan",
                 "Omnivore", "Omnivore", "Omnivore"),
  row.names  = c("SRR8146968", "SRR8146978", "SRR8146982",
                 "SRR8146938", "SRR8146935", "SRR8146936")
)

# Build phyloseq object
OTU  <- otu_table(as.matrix(otu_counts), taxa_are_rows = TRUE)
META <- sample_data(metadata)
ps   <- phyloseq(OTU, META)


# 2. TAXONOMIC ABUNDANCE — PHYLUM LEVEL  -----------------
# Bracken outputs species-level data. We map species to phylum using a
# lookup based on known gut microbiome taxonomy, then aggregate.

# Convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps_melt <- psmelt(ps_rel)

# Map species to phylum using genus-level keywords
assign_phylum <- function(species) {
  s <- tolower(species)
  case_when(
    grepl("bacteroides|parabacteroides|alistipes|phocaeicola|prevotella|paraprevotella|segatella|odoribacter|barnesiella", s) ~ "Bacteroidota",
    grepl("faecalibacterium|ruminococcus|lachnospira|coprococcus|agathobacter|butyrivibrio|roseburia|blautia|eubacterium|dorea|subdoligranulum|flavonifractor|oscillibacter|lawsonibacter|dysmobacter|dysosmobacter|wujia|simiaoa", s) ~ "Firmicutes",
    grepl("bifidobacterium|collinsella|slackia", s) ~ "Actinobacteriota",
    grepl("akkermansia", s) ~ "Verrucomicrobiota",
    grepl("sutterella|bilophila", s) ~ "Proteobacteria",
    grepl("fusobacterium", s) ~ "Fusobacteriota",
    grepl("escherichia|citrobacter|enterobacter|klebsiella", s) ~ "Proteobacteria",
    TRUE ~ "Other"
  )
}

ps_melt <- ps_melt %>%
  mutate(Phylum = assign_phylum(OTU)) %>%
  group_by(Sample, Phylum) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  left_join(metadata %>% tibble::rownames_to_column("Sample"), by = "Sample") %>%
  # Clean sample labels
  mutate(Sample = recode(Sample,
                         "SRR8146968" = "Vegan_1",
                         "SRR8146978" = "Vegan_2",
                         "SRR8146982" = "Vegan_3",
                         "SRR8146938" = "Omnivore_1",
                         "SRR8146935" = "Omnivore_2",
                         "SRR8146936" = "Omnivore_3"
  ))

# Color palette for phyla
phylum_colors <- c(
  "Bacteroidota"       = "#4E79A7",
  "Firmicutes"         = "#F28E2B",
  "Actinobacteriota"   = "#59A14F",
  "Verrucomicrobiota"  = "#B07AA1",
  "Proteobacteria"     = "#E15759",
  "Fusobacteriota"     = "#76B7B2",
  "Other"              = "#BAB0AC"
)

p_abundance <- ggplot(ps_melt, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ diet_group, scales = "free_x") +
  scale_fill_manual(values = phylum_colors) +
  labs(title = "Taxonomic Abundance by Sample (Phylum Level)",
       x = "Sample", y = "Relative Abundance",
       fill = "Phylum") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave("figures/taxonomic_abundance.png", p_abundance, width = 10, height = 7, dpi = 300)
cat("Saved: figures/taxonomic_abundance.pdf\n")


# 3. ALPHA DIVERSITY  -----------------
# Calculate alpha diversity metrics
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Sample    <- rownames(alpha_div)
alpha_div$diet_group <- metadata[alpha_div$Sample, "diet_group"]
alpha_div$Sample <- recode(alpha_div$Sample,
                           "SRR8146968" = "Vegan_1", "SRR8146978" = "Vegan_2", "SRR8146982" = "Vegan_3",
                           "SRR8146938" = "Omnivore_1", "SRR8146935" = "Omnivore_2", "SRR8146936" = "Omnivore_3"
)

# Print summary
cat("\nAlpha Diversity Summary:\n")
print(alpha_div)

# Wilcoxon test between groups
cat("\nWilcoxon test - Shannon diversity:\n")
print(wilcox.test(Shannon ~ diet_group, data = alpha_div))

cat("\nWilcoxon test - Observed richness:\n")
print(wilcox.test(Observed ~ diet_group, data = alpha_div))

# Plot alpha diversity
alpha_long <- tidyr::pivot_longer(alpha_div,
                                  cols = c("Observed", "Shannon", "Simpson"),
                                  names_to = "Metric", values_to = "Value")

p_alpha <- ggplot(alpha_long, aes(x = diet_group, y = Value, fill = diet_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  facet_wrap(~ Metric, scales = "free_y") +
  scale_fill_manual(values = c("Vegan" = "#2ecc71", "Omnivore" = "#e74c3c")) +
  labs(title = "Alpha Diversity by Dietary Group",
       x = "Diet Group", y = "Diversity Value", fill = "Diet Group") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("figures/alpha_diversity.png", p_alpha, width = 8, height = 5, dpi = 300)
cat("Saved: figures/alpha_diversity.pdf\n")


# 4. BETA DIVERSITY  -----------------
# Calculate Bray-Curtis dissimilarity on relative abundance
otu_rel <- sweep(otu_counts, 2, colSums(otu_counts), "/")
bray_dist <- vegdist(t(otu_rel), method = "bray")

# PERMANOVA test
cat("\nPERMANOVA (Bray-Curtis, diet_group):\n")
permanova <- adonis2(bray_dist ~ diet_group, data = metadata, permutations = 999)
print(permanova)

# PCoA ordination
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  PC1        = pcoa$points[, 1],
  PC2        = pcoa$points[, 2],
  Sample     = rownames(pcoa$points),
  diet_group = metadata[rownames(pcoa$points), "diet_group"]
) %>% mutate(Sample = recode(Sample,
                             "SRR8146968" = "Vegan_1", "SRR8146978" = "Vegan_2", "SRR8146982" = "Vegan_3",
                             "SRR8146938" = "Omnivore_1", "SRR8146935" = "Omnivore_2", "SRR8146936" = "Omnivore_3"
))

# Percent variance explained
var_explained <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = diet_group, label = Sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3) +
  scale_color_manual(values = c("Vegan" = "#2ecc71", "Omnivore" = "#e74c3c")) +
  labs(title = "Beta Diversity — Bray-Curtis PCoA",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)"),
       color = "Diet Group") +
  theme_bw() +
  annotate("text", x = Inf, y = -Inf,
           label = paste0("PERMANOVA p = ", permanova$`Pr(>F)`[1]),
           hjust = 1.1, vjust = -0.5, size = 3.5)

ggsave("figures/beta_diversity_pcoa.png", p_beta, width = 7, height = 6, dpi = 300)
cat("Saved: figures/beta_diversity_pcoa.pdf\n")

# Save OTU counts for DESeq2
saveRDS(list(otu_counts = otu_counts, metadata = metadata),
        file = "results/bracken/data_for_deseq2.rds")
cat("Saved count data for DESeq2: results/bracken/data_for_deseq2.rds\n")
