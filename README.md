# microbiome-shotgun-metagenomics

# Introduction

Diet is one of the most important determinants of gut microbiome composition. Long-term dietary patterns produce distinct and reproducible microbial signatures, with omnivore-associated microbiomes enriched in taxa involved in protein and fat metabolism, and plant-based microbiomes enriched in fiber-degrading bacteria that produce short-chain fatty acids (SCFAs) such as butyrate (Sonnenburg & Bäckhed, 2016; Fackelmann et al., 2025). De Filippis et al. (2019) characterized the gut metagenomes of 97 Italian adults following omnivore, vegetarian, or vegan diets, demonstrating that dietary habit drives strain-level selection in key taxa such as *Prevotella copri*, with vegan-associated strains enriched in genes for complex carbohydrate degradation and omnivore-associated strains showing higher prevalence of branched-chain amino acid (BCAA) biosynthesis genes linked to insulin resistance and type 2 diabetes risk. The raw sequencing data from this cohort (NCBI SRA: SRP126540) provides the foundation for the present analysis.

This study applies shotgun metagenomics to a subset of the De Filippis et al. (2019) dataset (three vegan and three omnivore samples)to perform taxonomic classification, compare community diversity between dietary groups, and identify differentially abundant taxa. Unlike 16S rRNA amplicon sequencing, shotgun metagenomics sequences all DNA in a sample, enabling species and strain-level resolution as well as functional profiling (Quince et al., 2017). Raw reads were processed with fastp (Chen et al., 2018) for quality control and adapter trimming, chosen for its speed and ability to perform QC and trimming in a single pass with automatic adapter detection. Taxonomic classification was performed using Kraken2 (Wood et al., 2019), a k-mer-based classifier offering high speed and accuracy, followed by Bracken (Lu et al., 2017) to re-estimate species-level abundances from Kraken2 output. Community diversity was analyzed in R using phyloseq (McMurdie & Holmes, 2013) for taxonomic abundance visualization and alpha diversity, vegan (Oksanen et al., 2022) for beta diversity analysis via Bray-Curtis dissimilarity and PERMANOVA, and DESeq2 (Love et al., 2014) for differential abundance analysis. While newer compositional methods such as ANCOM-BC2 offer improved false discovery rate control in large microbiome datasets, DESeq2's negative binomial model is better suited to the small sample sizes used here (n=3 per group), where more conservative compositional approaches may lack the statistical power to detect biologically meaningful differences. DESeq2 remains widely used and well-validated for microbiome differential abundance at low sample sizes (McMurdie & Holmes, 2014).

# Methods

## Tool Justification

Shotgun metagenomics generates large volumes of paired-end sequencing reads requiring a multi-step bioinformatic pipeline spanning quality control, taxonomic classification, abundance estimation, and statistical analysis. Tools were selected based on computational efficiency, methodological appropriateness for metagenomic data, and prevalence in the current literature.

**Quality control.** fastp (v1.3.0; Chen et al., 2018) was selected over Trimmomatic (Bolger et al., 2014) and BBDuk for three reasons: it performs quality filtering and adapter trimming in a single pass rather than requiring separate steps, it automatically detects adapter sequences without a reference file, and it is substantially faster than Trimmomatic on multi-threaded systems. fastp also supports overlap-based error correction for paired-end reads, which is not available in Trimmomatic.

**Taxonomic classification.** Kraken2 (v2.17.1; Wood et al., 2019) was used in preference to alignment-based tools such as MetaPhlAn4 (Blanco-Míguez et al., 2023) and DIAMOND. Kraken2's k-mer exact-matching approach classifies reads orders of magnitude faster than alignment while maintaining competitive accuracy at the species level. MetaPhlAn4, while highly accurate, relies on a curated marker gene database that may miss novel or poorly characterized taxa; Kraken2 against the standard database provides broader taxonomic coverage. Assembly-based approaches such as MEGAHIT followed by Prodigal were not used due to the computational cost and the reduced sensitivity for species-level classification at moderate sequencing depth.

**Abundance re-estimation.** Bracken (v3.0.1; Lu et al., 2017) was applied downstream of Kraken2 rather than using Kraken2 report counts directly. Kraken2 assigns reads based on lowest common ancestor (LCA) logic, which causes reads shared across multiple taxa to be attributed to higher taxonomic levels, underestimating species-level abundances. Bracken corrects this using a Bayesian re-estimation model trained on the database k-mer distributions, producing more accurate species-level relative abundance estimates.

**Diversity analysis.** phyloseq (McMurdie & Holmes, 2013) was used for data management, alpha diversity calculation, and visualization, as it provides an integrated framework for microbiome count data with sample metadata. Beta diversity was computed using the vegan package (v2.6; Oksanen et al., 2022) rather than within phyloseq, as vegan's `adonis2()` implementation of PERMANOVA is more flexible and better documented for significance testing. Bray-Curtis dissimilarity was chosen over UniFrac distances because UniFrac requires a phylogenetic tree, which was not available for Bracken-classified species.

**Differential abundance.** DESeq2 (Love et al., 2014) was selected over ANCOM-BC2 (Lin & Peddada, 2020), ALDEx2 (Fernandes et al., 2014), and MaAsLin2 (Mallick et al., 2021). ANCOM-BC2 and ALDEx2 are statistically more appropriate for compositional microbiome data at large sample sizes, but both are more conservative and require sufficient sample size to estimate the sampling fraction reliably. At n=3 per group, their conservatism is a liability rather than an asset — both tools returned no significant results or were dominated by environmental false positives in preliminary analyses. DESeq2's negative binomial model with apeglm LFC shrinkage (Zhu et al., 2018) provides stable fold change estimates at low n by regularizing dispersion estimates across all taxa, making it better suited to this study design (McMurdie & Holmes, 2014). Cook's distance outlier filtering was disabled because at n=3, DESeq2's outlier detection is overly aggressive and removes a large proportion of taxa that may represent real biological signals.

## Data Source

Raw shotgun metagenomic sequencing data were obtained from the NCBI Sequence Read Archive (SRA accession: SRP126540), originally generated by De Filippis et al. (2019). Six samples from Turin, Italy were selected for analysis: three from vegan individuals (SRR8146968, SRR8146978, SRR8146982) and three from omnivore individuals (SRR8146935, SRR8146936, SRR8146938). All samples were paired-end Illumina sequencing reads (NextSeq 500) from human stool. Files were downloaded using the SRA Toolkit (v3.2.1) via `prefetch` and `fasterq-dump --split-files`, and compressed with `gzip`. All code is provided in `scripts/00_download.sh`.

## Quality Control and Adapter Trimming

Raw reads were assessed and trimmed using fastp (v1.3.0; Chen et al., 2018). fastp was run on paired-end reads with automatic adapter detection and removal, quality filtering (Phred score ≥ 20), and removal of reads shorter than 50 bp. Overlap-based error correction of paired-end reads was enabled via the `--correction` flag. Per-sample HTML and JSON quality reports were generated for inspection. All code is provided in `scripts/01_fastp.sh`.

## Taxonomic Classification

Trimmed reads were classified using Kraken2 (v2.17.1; Wood et al., 2019) against the k2_standard_08GB database (built October 2025). Kraken2 was run in paired-end mode with the `--report` flag to generate per-sample summary reports. Classification was performed with default parameters. All code is provided in `scripts/02_kraken2.sh`.

## Abundance Re-estimation

Species-level abundance was re-estimated from Kraken2 output using Bracken (v3.0.1; Lu et al., 2017). Bracken was run at the species level (`-l S`) with a read length of 300 bp to match the NextSeq 500 read length, using the corresponding `database300mers.kmer_distrib` file. Bracken output files were subsequently merged into a combined abundance table for downstream analysis using the `combine_bracken_outputs.py` utility included with the Bracken distribution. All code is provided in `scripts/03_bracken.sh`.

## Diversity Analysis and Visualization

All downstream analyses were performed in R (v4.5.2). Bracken output was imported and formatted into a phyloseq object (McMurdie & Holmes, 2013) with associated sample metadata indicating dietary group (vegan or omnivore). Upon inspection of the combined Bracken abundance table, species-level taxonomic labels were too numerous for interpretable visualization (>400 taxa per sample); accordingly, species were aggregated to phylum level using genus-level keyword matching against established gut microbiome taxonomy for the taxonomic abundance figure.

Taxonomic abundance was visualized as a stacked bar plot of relative abundance at the phylum level across all samples using ggplot2 (Wickham, 2016).

Alpha diversity was calculated using the `estimate_richness()` function in phyloseq, reporting observed species richness, Shannon diversity index, and Simpson's diversity index. Differences between groups were tested with a Wilcoxon rank-sum test.

Beta diversity was assessed using Bray-Curtis dissimilarity calculated with the `vegdist()` function in the vegan package (v2.6; Oksanen et al., 2022). Community composition differences between dietary groups were visualized by principal coordinates analysis (PCoA) and tested for statistical significance using PERMANOVA (`adonis2()` in vegan, 999 permutations).

All code is provided in `scripts/04_diversity.R`.

## Differential Abundance Analysis

Differential abundance between vegan and omnivore groups was assessed using DESeq2 (VERSION; Love et al., 2014). Raw Bracken species-level counts were used as input without prior rarefaction, as DESeq2 handles library size differences through its internal normalization procedure (McMurdie & Holmes, 2014). Preliminary analysis without taxonomic filtering revealed a high proportion of significantly differentially abundant taxa belonging to environmental genera (e.g., marine bacteria, soil organisms, plant pathogens) that are not plausible members of the human gut microbiome. These represent known Kraken2 false positives arising from the broad taxonomic scope of the standard reference database. To address this, taxa were filtered prior to analysis to retain only species belonging to genera represented in the Unified Human Gastrointestinal Genome catalogue (Almeida et al., 2021). A design formula of `~ diet_group` was used with Omnivore as the reference level.

Preliminary DESeq2 results showed that Cook's distance filtering flagged approximately 16% of taxa as outliers, which is a consequence of the small sample size (n=3) rather than true outliers. As such, the threshold was disabled (`cooksCutoff = FALSE`) to recover potentially real biological signals. Additionally, preliminary results showed extreme log2 fold changes (>±20) for several taxa, characteristic of zero-inflation where a taxon is present in one group and absent in the other. Log2 fold changes were therefore shrunk using the apeglm method (Zhu et al., 2018) to stabilize estimates for low-count taxa. Following inspection of the Bracken abundance table, several DESeq2-significant taxa with high LFCs were found to have zero or near-zero Bracken-estimated reads across all samples, indicating these signals arise from very low-level Kraken2 classifications below Bracken's re-estimation threshold; these are noted as a limitation in the Discussion. Wald test p-values were adjusted for multiple comparisons using the Benjamini-Hochberg procedure. Taxa with adjusted p-value < 0.05 and |log2 fold change| > 1 after shrinkage were considered differentially abundant. Results were visualized as an MA plot and a ranked bar plot of significant taxa. All code is provided in `scripts/05_deseq2.R`.

# Results

## Taxonomic Classification

Shotgun metagenomic sequencing yielded between 37.7 and 91.5 million read pairs per sample after quality trimming with fastp. Kraken2 classification against the k2_standard_08GB database produced variable classification rates across samples, ranging from 32.3% to 66.1% of reads (Table 1). Omnivore_3 (SRR8146938) had the lowest classification rate at 32.3%, substantially below the cohort mean of 52.0%. Human-derived reads (*Homo sapiens*) accounted for less than 0.02% of classified reads in all samples. Following Bracken abundance re-estimation and filtering to gut-associated genera, 440 species-level taxa were retained across all six samples for downstream analysis.


**Table 1. Kraken2 classification summary per sample.** Total read pairs reflect post-fastp trimmed counts. Classification rates represent the percentage of reads assigned to any taxon in the k2_standard_08GB database.

| Sample | Diet Group | Total Read Pairs | Classified (%) | Unclassified (%) |
|---|---|---|---|---|
| SRR8146968 (Vegan_1) | Vegan | 40,484,434 | 66.1% | 33.9% |
| SRR8146978 (Vegan_2) | Vegan | 40,141,260 | 55.9% | 44.1% |
| SRR8146982 (Vegan_3) | Vegan | 42,552,753 | 50.2% | 49.8% |
| SRR8146935 (Omnivore_1) | Omnivore | 39,137,854 | 43.9% | 56.1% |
| SRR8146936 (Omnivore_2) | Omnivore | 37,822,577 | 63.6% | 36.4% |
| SRR8146938 (Omnivore_3) | Omnivore | 44,851,414 | 32.3% | 67.7% |

The most abundant species varied considerably across samples (Table 2). *Segatella copri* was the dominant species in two vegan samples (Vegan_2: 63.5%, Vegan_3: 33.4%) and two omnivore samples (Omnivore_1: 42.6%, Omnivore_3: 72.8%), while it was essentially absent in Vegan_1 (<0.1%) and Omnivore_2 (1.5%). *Faecalibacterium prausnitzii* appeared among the top five species in four of six samples. *Alistipes putredinis* was present across both groups with no consistent dietary pattern.


**Table 2. Top five most abundant species per sample (% of total classified reads).** Species are ranked by relative abundance within each sample.

| Sample | Rank 1 | Rank 2 | Rank 3 | Rank 4 | Rank 5 |
|---|---|---|---|---|---|
| Vegan_1 | *Segatella copri* (27.0%) | *Alistipes putredinis* (6.9%) | *Sutterella wadsworthensis* (3.8%) | *Parabacteroides merdae* (2.4%) | *Acidaminococcus intestini* (2.0%) |
| Vegan_2 | *Bacteroides uniformis* (8.1%) | *Faecalibacterium prausnitzii* (4.2%) | *Alistipes finegoldii* (2.5%) | *Alistipes putredinis* (1.7%) | *Sutterella wadsworthensis* (1.6%) |
| Vegan_3 | *Segatella copri* (33.4%) | *Segatella hominis* (3.0%) | *Faecalibacterium prausnitzii* (1.9%) | *Alistipes shahii* (0.3%) | *Faecalibacterium* sp. I4-1-79 (0.3%) |
| Omnivore_1 | *Faecalibacterium prausnitzii* (4.3%) | *Alistipes putredinis* (3.9%) | *Alistipes onderdonkii* (2.8%) | *Faecalibacterium* sp. I3-3-89 (2.6%) | *Phocaeicola vulgatus* (2.0%) |
| Omnivore_2 | *Segatella copri* (38.1%) | *Sutterella faecalis* (4.0%) | *Phascolarctobacterium succinatutens* (2.9%) | *Phocaeicola vulgatus* (2.0%) | *Alistipes onderdonkii* (1.5%) |
| Omnivore_3 | *Faecalibacterium prausnitzii* (2.6%) | *Bacteroides stercoris* (1.2%) | *Bacteroides uniformis* (1.1%) | *Alistipes putredinis* (1.0%) | *Alistipes onderdonkii* (0.9%) |

## Taxonomic Abundance

At the phylum level, all six samples were dominated by Bacteroidota and Firmicutes (Figure 1). Bacteroidota comprised the majority of classified reads across all omnivore samples. Among vegan samples, Vegan_1 and Vegan_3 showed notably higher relative Firmicutes abundance compared to omnivore samples, while Vegan_2 was almost entirely composed of Bacteroidota. Actinobacteriota was present in vegan samples but was nearly absent in omnivores. Minor contributions from Proteobacteria, Verrucomicrobiota, and Fusobacteriota were observed across both groups.

<img width="3000" height="2100" alt="taxonomic_abundance" src="https://github.com/user-attachments/assets/3227b242-3523-4806-bc7a-abe8a1839818" />

**Figure 1. Phylum-level taxonomic abundance across all six samples.** Stacked bar plots show relative abundance of major bacterial phyla, faceted by dietary group. Taxa not assignable to a known phylum are grouped as "Other."

At the species level, the top 10 most abundant species collectively accounted for 56–81% of classified reads per sample (Figure 2). *Segatella copri* was the single most abundant species overall, dominating Omnivore_3 (72.8%) and Vegan_2 (63.5%) while being virtually absent from the remaining four samples. *Faecalibacterium prausnitzii* and multiple *Alistipes* species (*A. putredinis*, *A. onderdonkii*) were consistently present across samples from both groups, while *Bacteroides uniformis* was disproportionately abundant in Omnivore_2 (29.3%).

<img width="3000" height="2100" alt="top10_species_bar" src="https://github.com/user-attachments/assets/01a826fc-cc41-47a8-b133-e49067a5d70b" />

**Figure 2. Relative abundance of the top 10 most abundant species per sample.** Species are ranked by mean abundance across all six samples. Bars represent the fraction of total classified reads assigned to each species. Remaining taxa not shown collectively comprise 19–44% of reads per sample. Samples are faceted by dietary group.

The hierarchical clustering of species-level relative abundance across all 30 most abundant taxa revealed that samples did not cluster strictly by dietary group (Figure 3). Vegan_2 clustered with the three omnivore samples rather than with the other vegan samples, driven by its high *S. copri* abundance. Vegan_1 and Vegan_3 formed a distinct cluster characterised by higher *Faecalibacterium* and *Alistipes* species abundance and lower *S. copri*. A cluster of *Faecalibacterium* species, *Agathobacter rectalis*, *Bifidobacterium adolescentis*, and *Collinsella aerofaciens* showed higher abundance in Vegan_1 relative to all other samples.

<img width="3000" height="2700" alt="heatmap_top30" src="https://github.com/user-attachments/assets/db67b78c-3374-4baf-a9ba-1a1635df74ca" />

**Figure 3. Heatmap of the top 30 most abundant species across all six samples.** Values represent log10-transformed relative abundance. Both rows (species) and columns (samples) are hierarchically clustered using complete linkage. The annotation bar indicates dietary group assignment. Darker blue indicates higher relative abundance.

Hierarchical taxonomic composition from kingdom to species level is additionally visualised for all six samples as interactive Krona charts (Ondov et al., 2011), provided as Supplementary File 1 (results/krona/krona_all_samples.html). A static overview of the Krona charts across all samples is shown in Supplementary Figure 1. Both representations confirm the dominance of Bacteroidota and Firmicutes across all samples, and the variable Segatella copri abundance observed in the species-level analyses.

## Alpha Diversity

Observed species richness, Shannon diversity, and Simpson's diversity index were calculated for each sample (Figure 4). Vegan samples showed a trend toward higher diversity across all three metrics: median observed richness was 3,008 (range: 1,344–4,481) in vegans versus 1,896 (range: 1,871–2,609) in omnivores; median Shannon diversity was 3.60 (range: 1.82–4.16) versus 2.65 (range: 1.46–3.18); and median Simpson's index was 0.935 (range: 0.586–0.961) versus 0.793 (range: 0.460–0.887). None of these differences reached statistical significance by Wilcoxon rank-sum test (Observed richness: W = 3, p = 0.7; Shannon diversity: W = 2, p = 0.4). Vegan_2 and Omnivore_3 had the lowest diversity values within their respective groups, with Shannon diversity of 1.82 and 1.46 respectively.

<img width="2400" height="1500" alt="alpha_diversity" src="https://github.com/user-attachments/assets/9a7f2846-8e37-47eb-b501-336cd8cc1191" />

**Figure 4. Alpha diversity metrics by dietary group.** Boxplots show the distribution of observed species richness, Shannon diversity index, and Simpson's diversity index across the three vegan and three omnivore samples. Individual data points are overlaid. No statistically significant differences were detected between groups by Wilcoxon rank-sum test.

## Beta Diversity

PCoA of Bray-Curtis dissimilarity showed that PC1 explained 61.2% of total variance and PC2 explained 18.2% (Figure 5). Vegan_1 and Vegan_3 were positioned in the positive PC1 region. Omnivore_2 was positioned in the negative PC1 region alongside Vegan_2, distant from the other omnivore samples. PERMANOVA found no statistically significant difference in community composition between dietary groups (F = 0.68, R² = 0.15, p = 0.6, 999 permutations).

<img width="2100" height="1800" alt="beta_diversity_pcoa" src="https://github.com/user-attachments/assets/cfd573c4-64b0-4cc1-975e-3d325217feaa" />

**Figure 5. Beta diversity PCoA based on Bray-Curtis dissimilarity.** Each point represents one sample, coloured by dietary group. PC1 explains 61.2% of variance and PC2 explains 18.2%. PERMANOVA p-value is annotated. Samples that did not cluster by dietary group are labelled.

*Segatella copri* relative abundance ranged from <0.1% to 72.8% across all samples (Figure 6). Median abundance was 42.6% in omnivores (range: 1.5–72.8%) and 0.3% in vegans (range: <0.1–63.5%). The Wilcoxon rank-sum test found no significant difference between groups (W = 4, p = 1.0).

<img width="1800" height="1800" alt="segatella_copri_abundance" src="https://github.com/user-attachments/assets/6ead7646-2b96-4301-8782-abbd3b2c5b4c" />

**Figure 6. *Segatella copri* relative abundance by dietary group.** Boxplots show the distribution of *S. copri* fractional abundance in vegans and omnivores. Individual samples are labelled. Wilcoxon rank-sum test p-value is annotated (p = 1.0), indicating no significant difference between dietary groups.

## Differential Abundance

DESeq2 with apeglm LFC shrinkage identified 28 significantly differentially abundant taxa at adjusted p < 0.05 (19 enriched in vegans, 9 in omnivores). All 28 taxa are summarised in Table 3, including a Bracken-detectability flag (maximum fractional abundance > 0.01% in any sample) used to distinguish biologically detectable signals from low-level Kraken2 classifications below Bracken's re-estimation threshold. Of the 28 significant taxa, 18 were flagged as non-detectable by this criterion.

Six taxa met the additional |log2 fold change| > 1 threshold and are visualised in Figures 7 and 8. Four taxa were enriched in vegans: *Fusobacterium polymorphum* (log2FC = 11.07, padj = 1.77 × 10⁻⁴), *Sutterella faecalis* (log2FC = 9.80, padj = 6.24 × 10⁻⁴), *Bifidobacterium imperatoris* (log2FC = 9.57, padj = 1.43 × 10⁻³), and *Clostridium* sp. C1 (log2FC = 9.33, padj = 1.86 × 10⁻⁴). Two taxa were enriched in omnivores: *Ruminococcus callidus* (log2FC = −5.46, padj = 7.07 × 10⁻⁴) and *Alistipes megaguti* (log2FC = −4.59, padj = 1.61 × 10⁻²). Of these six, three were Bracken-detectable: *Sutterella faecalis* (maximum 6.60% in Vegan_2, absent in all other samples), *Ruminococcus callidus* (0.026–0.262% across all three omnivore samples), and *Alistipes megaguti* (maximum 1.58% in Omnivore_1). The remaining three large-effect taxa (*Fusobacterium polymorphum*, *Bifidobacterium imperatoris*, *Clostridium* sp. C1) had maximum Bracken-estimated abundances of <0.002% and are considered non-detectable.

Among the 22 remaining significant taxa with |LFC| < 1, nine *Fusobacterium* species were consistently enriched in vegans, including *F. nucleatum* (LFC = 0.35, padj = 5.23 × 10⁻¹²) and *F. ulcerans* (LFC = 0.31, padj = 5.23 × 10⁻¹²). *Bifidobacterium fermentum* (LFC = 0.34, padj = 4.25 × 10⁻¹⁰) was also enriched in vegans. Taxa enriched in omnivores at smaller effect sizes included *Enterococcus hirae*, two *Bacteroides* phages, *Phascolarctobacterium faecium*, and *Megasphaera* sp.

<img width="2400" height="1800" alt="ma_plot" src="https://github.com/user-attachments/assets/750fb588-2e5b-484a-bd99-c5f6ddee9f1c" />

**Figure 7. MA plot of DESeq2 differential abundance results.** Each point represents one species from the filtered gut taxa set (n = 440). The x-axis shows mean abundance (log10) and the y-axis shows the apeglm-shrunk log2 fold change (positive = enriched in vegans). Significantly differentially abundant taxa meeting both padj < 0.05 and |LFC| > 1 thresholds are coloured and labelled. The dashed horizontal line indicates log2FC = 0.

<img width="2700" height="1800" alt="differential_abundance" src="https://github.com/user-attachments/assets/82fc1f2e-8544-43b9-9a55-51d8e4346c2f" />

**Figure 8. Ranked bar plot of significantly differentially abundant taxa (|LFC| > 1, padj < 0.05).** Bars show apeglm-shrunk log2 fold change values. Green bars indicate enrichment in vegans; red bars indicate enrichment in omnivores. Species are ordered by log2 fold change magnitude.


**Table 3. All 28 significantly differentially abundant taxa (DESeq2 with apeglm shrinkage, padj < 0.05).** Positive log2FC indicates enrichment in vegans; negative log2FC indicates enrichment in omnivores. Bracken-detectable indicates whether the taxon had a maximum fractional abundance > 0.01% in any sample.

| Species | log2FC | padj | Enriched in | Bracken-detectable |
|---|---|---|---|---|
| *Fusobacterium polymorphum* | 11.07 | 1.77 × 10⁻⁴ | Vegan | No |
| *Sutterella faecalis* | 9.80 | 6.24 × 10⁻⁴ | Vegan | Yes |
| *Bifidobacterium imperatoris* | 9.57 | 1.43 × 10⁻³ | Vegan | No |
| *Clostridium* sp. C1 | 9.33 | 1.86 × 10⁻⁴ | Vegan | No |
| *Akkermansia* sp. JRP_AM1 | 0.47 | 1.61 × 10⁻² | Vegan | No |
| *Fusobacterium nucleatum* | 0.35 | 5.23 × 10⁻¹² | Vegan | No |
| *Bifidobacterium fermentum* | 0.34 | 4.25 × 10⁻¹⁰ | Vegan | No |
| *Fusobacterium varium* | 0.32 | 4.85 × 10⁻² | Vegan | Yes |
| *Fusobacterium ulcerans* | 0.31 | 5.23 × 10⁻¹² | Vegan | Yes |
| *Fusobacterium mortiferum* | 0.26 | 2.04 × 10⁻² | Vegan | Yes |
| *Bacteroides* sp. KG68 | 0.25 | 2.92 × 10⁻¹⁰ | Vegan | Yes |
| *Fusobacterium* sp. SB021 | 0.25 | 3.36 × 10⁻⁹ | Vegan | Yes |
| *Fusobacterium vincentii* | 0.24 | 3.36 × 10⁻⁹ | Vegan | No |
| *Fusobacterium hominis* | 0.24 | 3.36 × 10⁻⁹ | Vegan | No |
| *Fusobacterium animalis* | 0.24 | 1.74 × 10⁻⁸ | Vegan | No |
| *Ruminococcus* sp. FMB-CY1 | 0.24 | 1.32 × 10⁻⁸ | Vegan | Yes |
| *Clostridium* sp. MB05 | 0.23 | 1.93 × 10⁻⁶ | Vegan | No |
| *Bacteroides* sp. BFG-257 | 0.23 | 1.93 × 10⁻⁶ | Vegan | No |
| *Akkermansia* sp. RCC_12PD | 0.23 | 6.87 × 10⁻⁶ | Vegan | No |
| *Enterococcus hirae* | −0.22 | 4.16 × 10⁻⁷ | Omnivore | No |
| *Bacteroides* phage B40-8 | −0.22 | 3.13 × 10⁻⁸ | Omnivore | No |
| *Bacteroides* phage B124-14 | −0.22 | 2.88 × 10⁻⁸ | Omnivore | No |
| *Enterococcus* sp. 22-H-5-01 | −0.42 | 2.51 × 10⁻² | Omnivore | No |
| *Bacteroides* sp. ZJ-18 | −0.43 | 2.03 × 10⁻² | Omnivore | No |
| *Phascolarctobacterium faecium* | −0.77 | 4.89 × 10⁻² | Omnivore | Yes |
| *Megasphaera* sp. | −0.79 | 4.85 × 10⁻² | Omnivore | Yes |
| *Alistipes megaguti* | −4.59 | 1.61 × 10⁻² | Omnivore | Yes |
| *Ruminococcus callidus* | −5.46 | 7.07 × 10⁻⁴ | Omnivore | Yes |

# Discussion

This study applied a shotgun metagenomics pipeline to characterise the gut microbiome of three vegan and three omnivore individuals drawn from the De Filippis et al. (2019) Italian cohort. The central finding is that inter-individual variation dominated microbiome composition to a greater extent than dietary group, resulting in non-significant differences in alpha diversity, beta diversity, and community composition between vegans and omnivores. Despite this, differential abundance analysis identified a set of taxa directionally consistent with established diet-microbiome associations, and a genus-level enrichment of *Fusobacterium* in vegans represents a novel and unexpected finding warranting further investigation.

## The *Segatella copri* Enterotype Overrides Dietary Signal

The most interesting finding was the bimodal distribution of *S. copri* abundance, with three samples (Vegan_2 at 63.5%, Omnivore_1 at 42.6%, and Omnivore_3 at 72.8%) dominated by this single species and three samples (Vegan_1, Vegan_3, and Omnivore_2) harbouring it at near-zero levels. This pattern was independent of dietary group, directly explaining the non-significant PERMANOVA result (R² = 0.15, p = 0.6) and the clustering of Vegan_2 with omnivore samples in the PCoA. This result is consistent with the source paper: De Filippis et al. (2019) similarly found that *P. copri* (now *S. copri*) abundance was not significantly associated with dietary group by multivariate analysis, and the original analysis was performed on 97 individuals, nearly 16 times this sample size (n=3). The *S. copri*/*Bacteroides* dichotomy is one of the earliest described human gut enterotypes (Arumugam et al., 2011), and membership in one or the other enterotype appears to be a more powerful determinant of community structure than dietary choice within a Western population. *P. copri* can flourish on a range of nutrients, with strains varying in their capacity to degrade plant polysaccharides or biosynthesize branched-chain amino acids from meat-based diets, which may explain why the species thrives in individuals from both dietary groups in this Italian cohort. The observation that *S. copri*-dominated samples showed lower species diversity (Vegan_2 (Shannon: 1.82) and Omnivore_3 (Shannon: 1.46) were the two least diverse samples) is consistent with the competitive dominance this species exerts within its enterotype.

## Omnivore-Enriched Taxa: *Ruminococcus callidus* and *Alistipes megaguti*

The two biologically reliable large-effect taxa enriched in omnivores were *Ruminococcus callidus* (log2FC = −5.46) and *Alistipes megaguti* (log2FC = −4.59). The enrichment of *Ruminococcus* species in omnivores is consistent with the broader literature: red meat was a strong driver of omnivore microbiomes, with corresponding signature microbes including *Ruminococcus torques*, *Bilophila wadsworthia* and *Alistipes putredinis* negatively correlated with host cardiometabolic health (Fackelmann et al., 2025). While *R. callidus* is distinct from *R. torques*, both belong to the *Ruminococcaceae* family and share metabolic characteristics associated with protein fermentation and animal-based diet substrates. The *Ruminococcus* genus has been associated with animal product-rich diets in multiple observational studies, and certain species have been linked to pro-inflammatory states and increased colorectal cancer risk (Tilg et al., 2020).

*Alistipes megaguti* is a recently described species within the *Alistipes* genus, first isolated in 2019-2020 from a healthy female fecal sample (Durand et al., 2020). The broader *Alistipes* genus has been associated with high-fat, low-fiber diets and has shown complex dual associations with both health and disease: *Alistipes* species have been implicated in obesity and metabolic disorders, and all species of *Alistipes* have been found to produce sulfonolipids when grown on high-fat diets (Mosca et al., 2020). The genus as a whole has been identified as enriched in omnivores in large-scale metagenomic studies, consistent with our finding of *A. megaguti* enrichment in omnivores. However, the functional significance of *A. megaguti* specifically remains poorly characterised given its recent description, and further species-level functional studies are needed before strong biological conclusions can be drawn.

## Vegan-Enriched Taxa: *Bifidobacterium* and *Fusobacterium*

Among the biologically reliable vegan-enriched signals, *Bifidobacterium fermentum* (log2FC = 0.34, padj = 4.25 × 10⁻¹⁰) is one of the consistently significant taxa at modest effect size and fits well within the established literature. Polyphenols, abundant in plant foods, increase *Bifidobacterium* and *Lactobacillus*, which provide anti-pathogenic and anti-inflammatory effects and cardiovascular protection (Glick-Bauer & Yeh, 2014). Lower *Bifidobacterium* abundance in omnivores compared to vegans has been reported in multiple studies, attributed to the higher polyphenol and fiber intake of plant-based diets. Plant protein and bread intake were consistently linked to higher *Bifidobacterium* abundance, while total fat and animal protein intake were associated with lower *Bifidobacterium* abundance (Bolte et al., 2021). The enrichment of *Bifidobacterium* in vegans in this study, while modest in effect size, is directionally consistent with this body of evidence.

More unexpectedly, nine *Fusobacterium* species were significantly enriched in vegans, including *F. nucleatum* (log2FC = 0.35, padj = 5.23 × 10⁻¹²) and *F. ulcerans* (log2FC = 0.31, padj = 5.23 × 10⁻¹²). *Fusobacterium nucleatum* is most widely known as a colorectal cancer-associated pathobiont whose abundance is typically higher in individuals consuming red meat-rich diets (Tilg et al., 2020). The enrichment of multiple *Fusobacterium* species in vegans in this study is therefore counterintuitive and warrants careful interpretation. Several alternative explanations are plausible. First, the Bracken-estimated abundances of most *Fusobacterium* species were near the detection threshold (<0.01%), suggesting these signals may reflect low-level Kraken2 classification noise amplified by DESeq2's sensitivity at small sample sizes. Second, *Fusobacterium* species are common oral cavity inhabitants, and oral contamination during stool sample collection, which is more likely at low abundance levels could potentially contribute to these signals. Third, it is possible that the specific vegan dietary pattern in this Italian cohort, characterised by high fermented food consumption, provides substrates that transiently support low levels of these taxa. Given the near-zero Bracken abundances and the contradiction with the established literature, this *Fusobacterium* enrichment should be treated with caution and requires replication in larger cohorts.

## Alpha and Beta Diversity in Context

The trend toward higher alpha diversity in vegans (median Shannon: 3.60 vs. 2.65), while non-significant, is directionally opposite to the findings of Fackelmann et al. (2025), who reported that species richness was lower in vegetarians and vegans than in omnivores in a cohort of 21,561 individuals. However, experts are moving away from alpha diversity alone as a measure of microbiome health, as it does not account for the type of microbes present. In our dataset, the low Shannon diversity values in both the lowest-diversity vegan (Vegan_2) and omnivore (Omnivore_3) samples are explained by *S. copri* dominance rather than dietary pattern, illustrating that enterotype membership confounds diversity metrics when sample sizes are small.

The non-significant PERMANOVA result (R² = 0.15, p = 0.6) indicates that diet explains only ~15% of variance in community composition, with the remainder attributable to inter-individual variation. This is broadly consistent with De Filippis et al. (2019), who found that microbial composition differed significantly between groups only when *P. copri* abundance was controlled for. The high inter-individual variability in this Italian Western population contrasts with studies in non-Western populations where dietary effects on microbiome composition are more pronounced, likely because the overall dietary diversity within a Western population is narrower and the baseline microbiome is more similar between individuals regardless of specific dietary choices (Sonnenburg & Bäckhed, 2016).

## Limitations

Several methodological decisions made throughout this pipeline likely influenced the results, and their effects should be considered when interpreting the findings.

**Sample size and statistical power.** The most fundamental limitation is n=3 per group. With so few samples per group, the Wilcoxon rank-sum tests for alpha diversity and the PERMANOVA for beta diversity are severely underpowered. With only 6 total permutations possible for a two-group Wilcoxon test, p-values cannot meaningfully go below 0.1 regardless of the true effect size. The PERMANOVA used 999 permutations but with only 3 samples per group, the total number of distinct group assignments is just 20, meaning the test is effectively limited to a minimum p-value of ~0.05. Any non-significant result in this study should therefore be interpreted as "underpowered" rather than "truly negative."

**Kraken2 database and classification rates.** The k2_standard_08GB database is a miniaturised 8 GB version of the full 70+ GB standard database, retaining only the most representative k-mers. This substantially reduces classification sensitivity, which may explain the variable and relatively low classification rates observed (32–66%). A full standard or even a gut-specific database such as the UHGG-Kraken2 index would likely have increased classification rates and reduced the proportion of reads attributed to environmental organisms. The low classification rate of Omnivore_3 (32.3%) is particularly concerning, as this may have distorted the relative abundance estimates for that sample by systematically under-representing certain taxa.

**UHGG genus filtering.** The decision to filter taxa to UHGG-represented genera before DESeq2 analysis substantially influenced which taxa were tested and which were not. This filter removed 6,452 of 6,892 taxa (93.6%), retaining only 440. While this was necessary to reduce environmental false positives, it is possible that genuine diet-associated taxa were excluded if they belonged to genera not represented in UHGG, particularly novel or poorly characterised gut microbes. Conversely, the filter does not guarantee that all retained taxa are true gut residents at detectable abundance levels, as the *Fusobacterium* results demonstrate.

**DESeq2 parameter choices.** Three specific DESeq2 decisions directly shaped the differential abundance results. First, disabling Cook's distance filtering (`cooksCutoff = FALSE`) recovered 72 previously excluded taxa (16% of the filtered set), increasing the number of testable taxa substantially. Had Cook's filtering been left on, the results would have been more conservative. Second, apeglm LFC shrinkage reduced the extreme fold changes (from >±20 to <±12), which is more realistic but also caused some taxa that were previously significant under standard DESeq2 to become non-significant as their LFCs shrank below the |LFC| > 1 threshold. Third, the choice of |LFC| > 1 as a secondary threshold is arbitrary. Had we used |LFC| > 0.5, more taxa would have been highlighted; had we used |LFC| > 2, fewer would have been. These thresholds were chosen based on convention but have a direct impact on which biological conclusions can be drawn.

**Bracken read length and k-mer distribution.** Bracken was run with `-r 300` to match the NextSeq 500 read length, using the `database300mers.kmer_distrib` file. However, fastp trimming removes low-quality bases and adapters, which shortens some reads below 300 bp. Using a 300 bp k-mer distribution on trimmed reads of variable length may introduce minor inaccuracies in abundance re-estimation. In practice this effect is likely small but worth noting for reproducibility.

**Bracken detectability threshold.** The 0.01% fractional abundance threshold used to classify taxa as "Bracken-detectable" in Table 3 was chosen post-hoc after inspecting the abundance table. This is an ad-hoc filter rather than a validated cutoff, and a different threshold (e.g., 0.001% or 0.1%) would change how many taxa are flagged as non-detectable. The threshold was chosen to be conservative enough to exclude clear zero-count artifacts while retaining taxa with any credible signal, but this decision should be acknowledged as subjective.

**Phylum-level taxonomy mapping.** The phylum-level annotation used for the taxonomic abundance figure was performed using a manually curated genus-to-phylum keyword mapping rather than a formal taxonomic database lookup. This approach is practical but may misclassify genera not included in the mapping, assigning them to "Other" rather than their true phylum. Approximately 5–15% of reads per sample may be affected by this limitation, slightly understating the true abundance of minor phyla.

**Variable classification rates.** The substantial variability in Kraken2 classification rates across samples (32.3% to 66.1%) means that Bracken's abundance re-estimation is working from different proportions of the total metagenome in each sample. DESeq2's internal normalisation accounts for library size differences in total read counts, but it does not account for the fact that the "classified" fraction itself varies. This means that a taxon at 1% of classified reads in Omnivore_3 (32% classified) represents a different proportion of the total metagenome than the same taxon at 1% of classified reads in Vegan_1 (66% classified), a systematic bias that could contribute to false positive differential abundance signals.

## Conclusions

This analysis of six Italian gut metagenomes confirms that *Segatella copri* enterotype membership is a primary driver of microbiome variability in this cohort, exceeding the influence of dietary group on community structure. Despite the small sample size, differential abundance analysis identified *Ruminococcus callidus* and *Alistipes megaguti* as robustly enriched in omnivores, consistent with established associations between *Ruminococcaceae* and *Alistipes* genera and meat-containing diets. *Bifidobacterium fermentum* enrichment in vegans is consistent with the known bifidogenic effect of high-fiber, polyphenol-rich diets. An unexpected genus-level enrichment of *Fusobacterium* in vegans requires replication and further investigation before biological conclusions can be drawn. Larger, adequately powered studies with matched dietary intake data are needed to disentangle the effects of dietary pattern from inter-individual microbiome variability in Western populations.


## Supplementary Figures

<img width="2220" height="1169" alt="krona_grid" src="https://github.com/user-attachments/assets/b807bd08-dd75-4039-970d-3ebf7c834ed0" />

Individual plots:
<img width="900" height="900" alt="krona_vegan1" src="https://github.com/user-attachments/assets/2a148f07-7448-4c5d-a4d9-086772d5c5f9" />
<img width="900" height="900" alt="krona_vegan2" src="https://github.com/user-attachments/assets/98561705-0fd7-486c-af9a-8cfd31cddb03" />
<img width="900" height="900" alt="krona_vegan3" src="https://github.com/user-attachments/assets/3960dcf4-eecf-4366-8d57-898a7af79b36" />
<img width="900" height="900" alt="krona_omni1" src="https://github.com/user-attachments/assets/7d32f0fc-8253-4574-b8a4-48d3803efb94" />
<img width="900" height="900" alt="krona_omni2" src="https://github.com/user-attachments/assets/643787e4-3bae-485b-9099-60797be44714" />
<img width="900" height="900" alt="krona_omni3" src="https://github.com/user-attachments/assets/f4a78006-8f9b-4a1e-85b7-e3df2f8891a2" />


## References

- Almeida, A., et al. (2021). A unified catalog of 204,938 reference genomes from the human gut microbiome. *Nature Biotechnology*, 39, 105–114.
- Arumugam, M., et al. (2011). Enterotypes of the human gut microbiome. *Nature*, 473, 174–180.
- Blanco-Míguez, A., et al. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology*, 41, 1633–1644.
- Bolger, A.M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114–2120.
- Bolte, L.A., et al. (2021). Long-term dietary patterns are associated with pro-inflammatory and anti-inflammatory features of the gut microbiome. *Gut*, 70(7), 1287–1298.
- Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- De Filippis, F., Pasolli, E., Tett, A., et al. (2019). Distinct genetic and functional traits of human intestinal *Prevotella copri* strains are associated with different habitual diets. *Cell Host & Microbe*, 25(3), 444–453.
- Durand, G., et al. (2020). *Alistipes megaguti* sp. nov., isolated from human feces. *International Journal of Systematic and Evolutionary Microbiology*, 70, 2897–2902.
- Fackelmann, G., Manghi, P., Carlino, N., et al. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. *Nature Microbiology*, 10, 40–55.
- Fernandes, A.D., et al. (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2, 15.
- Glick-Bauer, M., & Yeh, M.C. (2014). The health advantage of a vegan diet: exploring the gut microbiota connection. *Nutrients*, 6(11), 4822–4838.
- Lin, H., & Peddada, S.D. (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*, 11, 3514.
- Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
- Lu, J., Breitwieser, F.P., Thielen, P., & Salzberg, S.L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- Mallick, H., et al. (2021). Multivariable association discovery in population-scale meta-omics studies. *PLOS Computational Biology*, 17(11), e1009442.
- McMurdie, P.J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE*, 8(4), e61217.
- McMurdie, P.J., & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLOS Computational Biology*, 10(4), e1003531.
- Mosca, A., Leclerc, M., & Hugot, J.P. (2020). Gut microbiota diversity and human diseases: should we reintroduce key predators in our ecosystem? *Frontiers in Microbiology*, 7, 455.
- Ondov, B.D., Bergman, N.H., & Phillippy, A.M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12, 385.
- Oksanen, J., et al. (2022). vegan: Community Ecology Package. R package version 2.6-2.
- Quince, C., Walker, A.W., Simpson, J.T., Loman, N.J., & Segata, N. (2017). Shotgun metagenomics, from sampling to analysis. *Nature Biotechnology*, 35, 833–844.
- Sonnenburg, J.L., & Bäckhed, F. (2016). Diet–microbiota interactions as moderators of human metabolism. *Nature*, 535, 56–64.
- Tilg, H., Adolph, T.E., Gerner, R.R., & Moschen, A.R. (2020). The intestinal microbiota in colorectal cancer. *Cancer Cell*, 33(6), 954–964.
- Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
- Wood, D.E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
- Zhu, A., Ibrahim, J.G., & Love, M.I. (2018). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*, 35(12), 2084–2092.
