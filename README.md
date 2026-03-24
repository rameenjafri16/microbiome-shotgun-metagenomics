# microbiome-shotgun-metagenomics

# Introduction

Diet is one of the most important determinants of gut microbiome composition. Long-term dietary patterns produce distinct and reproducible microbial signatures, with omnivore-associated microbiomes enriched in taxa involved in protein and fat metabolism, and plant-based microbiomes enriched in fiber-degrading bacteria that produce short-chain fatty acids (SCFAs) such as butyrate (Sonnenburg & Bäckhed, 2016; Fackelmann et al., 2025). De Filippis et al. (2019) characterized the gut metagenomes of 97 Italian adults following omnivore, vegetarian, or vegan diets, demonstrating that dietary habit drives strain-level selection in key taxa such as *Prevotella copri*, with vegan-associated strains enriched in genes for complex carbohydrate degradation and omnivore-associated strains showing higher prevalence of branched-chain amino acid (BCAA) biosynthesis genes linked to insulin resistance and type 2 diabetes risk. The raw sequencing data from this cohort (NCBI SRA: SRP126540) provides the foundation for the present analysis.

This study applies shotgun metagenomics to a subset of the De Filippis et al. (2019) dataset (three vegan and three omnivore samples)to perform taxonomic classification, compare community diversity between dietary groups, and identify differentially abundant taxa. Unlike 16S rRNA amplicon sequencing, shotgun metagenomics sequences all DNA in a sample, enabling species- and strain-level resolution as well as functional profiling (Quince et al., 2017). Raw reads were processed with fastp (Chen et al., 2018) for quality control and adapter trimming, chosen for its speed and ability to perform QC and trimming in a single pass with automatic adapter detection. Taxonomic classification was performed using Kraken2 (Wood et al., 2019), a k-mer-based classifier offering high speed and accuracy, followed by Bracken (Lu et al., 2017) to re-estimate species-level abundances from Kraken2 output. Community diversity was analyzed in R using phyloseq (McMurdie & Holmes, 2013) for taxonomic abundance visualization and alpha diversity, vegan (Oksanen et al., 2022) for beta diversity analysis via Bray-Curtis dissimilarity and PERMANOVA, and DESeq2 (Love et al., 2014) for differential abundance analysis. While newer compositional methods such as ANCOM-BC2 offer improved false discovery rate control in large microbiome datasets, DESeq2's negative binomial model is better suited to the small sample sizes used here (n=3 per group), where more conservative compositional approaches may lack the statistical power to detect biologically meaningful differences. DESeq2 remains widely used and well-validated for microbiome differential abundance at low sample sizes (McMurdie & Holmes, 2014).

# Methods

## Tool Justification

Shotgun metagenomics generates large volumes of paired-end sequencing reads requiring a multi-step bioinformatic pipeline spanning quality control, taxonomic classification, abundance estimation, and statistical analysis. Tools were selected based on computational efficiency, methodological appropriateness for metagenomic data, and prevalence in the current literature.

**Quality control.** fastp (v1.3.0; Chen et al., 2018) was selected over Trimmomatic (Bolger et al., 2014) and BBDuk for three reasons: it performs quality filtering and adapter trimming in a single pass rather than requiring separate steps, it automatically detects adapter sequences without a reference file, and it is substantially faster than Trimmomatic on multi-threaded systems. fastp also supports overlap-based error correction for paired-end reads, which is not available in Trimmomatic.

**Taxonomic classification.** Kraken2 (v2.17.1; Wood et al., 2019) was used in preference to alignment-based tools such as MetaPhlAn4 (Blanco-Míguez et al., 2023) and DIAMOND. Kraken2's k-mer exact-matching approach classifies reads orders of magnitude faster than alignment while maintaining competitive accuracy at the species level. MetaPhlAn4, while highly accurate, relies on a curated marker gene database that may miss novel or poorly characterized taxa; Kraken2 against the standard database provides broader taxonomic coverage. Assembly-based approaches such as MEGAHIT followed by Prodigal were not used due to the computational cost and the reduced sensitivity for species-level classification at moderate sequencing depth.

**Abundance re-estimation.** Bracken (v3.0.1; Lu et al., 2017) was applied downstream of Kraken2 rather than using Kraken2 report counts directly. Kraken2 assigns reads based on lowest common ancestor (LCA) logic, which causes reads shared across multiple taxa to be attributed to higher taxonomic levels, underestimating species-level abundances. Bracken corrects this using a Bayesian re-estimation model trained on the database k-mer distributions, producing more accurate species-level relative abundance estimates.

**Diversity analysis.** phyloseq (McMurdie & Holmes, 2013) was used for data management, alpha diversity calculation, and visualization, as it provides an integrated framework for microbiome count data with sample metadata. Beta diversity was computed using the vegan package (v2.6; Oksanen et al., 2022) rather than within phyloseq, as vegan's `adonis2()` implementation of PERMANOVA is more flexible and better documented for significance testing. Bray-Curtis dissimilarity was chosen over UniFrac distances because UniFrac requires a phylogenetic tree, which was not available for Bracken-classified species.

**Differential abundance.** DESeq2 (Love et al., 2014) was selected over ANCOM-BC2 (Lin & Peddada, 2020), ALDEx2 (Fernandes et al., 2014), and MaAsLin2 (Mallick et al., 2021). ANCOM-BC2 and ALDEx2 are statistically more appropriate for compositional microbiome data at large sample sizes, but both are more conservative and require sufficient sample size to estimate the sampling fraction reliably. At n=3 per group, their conservatism is a liability rather than an asset. Both tools returned no significant results or were dominated by environmental false positives in preliminary analyses. DESeq2's negative binomial model with apeglm LFC shrinkage (Zhu et al., 2018) provides stable fold change estimates at low n by regularizing dispersion estimates across all taxa, making it better suited to this study design (McMurdie & Holmes, 2014). Cook's distance outlier filtering was disabled because at n=3, DESeq2's outlier detection is overly aggressive and removes a large proportion of taxa that may represent real biological signals.

## Data Source

Raw shotgun metagenomic sequencing data were obtained from the NCBI Sequence Read Archive (SRA accession: SRP126540), originally generated by De Filippis et al. (2019). Six samples from Turin, Italy were selected for analysis: three from vegan individuals (SRR8146968, SRR8146978, SRR8146982) and three from omnivore individuals (SRR8146935, SRR8146936, SRR8146938). All samples were paired-end Illumina sequencing reads (NextSeq 500) from human stool. Files were downloaded using the SRA Toolkit (v3.2.1) via `prefetch` and `fasterq-dump --split-files`, and compressed with `gzip`. All code is provided in `scripts/00_download.sh`.

## Quality Control and Adapter Trimming

Raw reads were assessed and trimmed using fastp (v1.3.0; Chen et al., 2018). fastp was run on paired-end reads with automatic adapter detection and removal, quality filtering (Phred score ≥ 20), and removal of reads shorter than 50 bp. Overlap-based error correction of paired-end reads was enabled via the `--correction` flag. Per-sample HTML and JSON quality reports were generated for inspection. All code is provided in `scripts/01_fastp.sh`.

## Taxonomic Classification

Trimmed reads were classified using Kraken2 (v2.17.1; Wood et al., 2019) against the k2_standard_08GB database (built October 2025). Kraken2 was run in paired-end mode with the `--report` flag to generate per-sample summary reports. Classification was performed with default parameters. All code is provided in `scripts/02_kraken2.sh`.

## Abundance Re-estimation

Species-level abundance was re-estimated from Kraken2 output using Bracken (v3.0.1; Lu et al., 2017). Bracken was run at the species level (`-l S`) with a read length of 300 bp to match the NextSeq 500 read length, using the corresponding `database300mers.kmer_distrib` file. Bracken output files were subsequently merged into a combined abundance table for downstream analysis using the `combine_bracken_outputs.py` utility included with the Bracken distribution. All code is provided in `scripts/03_bracken.sh`.

## Diversity Analysis and Visualization

All downstream analyses were performed in R (v4.x). Bracken output was imported and formatted into a phyloseq object (McMurdie & Holmes, 2013) with associated sample metadata indicating dietary group (vegan or omnivore). Species names were mapped to phylum-level taxonomy using genus-level keyword matching against established gut microbiome taxonomy for visualization purposes.

Taxonomic abundance was visualized as a stacked bar plot of relative abundance at the phylum level across all samples using ggplot2 (Wickham, 2016).

Alpha diversity was calculated using the `estimate_richness()` function in phyloseq, reporting observed species richness, Shannon diversity index, and Simpson's diversity index. Differences between groups were tested with a Wilcoxon rank-sum test.

Beta diversity was assessed using Bray-Curtis dissimilarity calculated with the `vegdist()` function in the vegan package (v2.6; Oksanen et al., 2022). Community composition differences between dietary groups were visualized by principal coordinates analysis (PCoA) and tested for statistical significance using PERMANOVA (`adonis2()` in vegan, 999 permutations).

All code is provided in `scripts/04_diversity.R`.

## Differential Abundance Analysis

Differential abundance between vegan and omnivore groups was assessed using DESeq2 (VERSION; Love et al., 2014). Raw Bracken species-level counts were used as input without prior rarefaction, as DESeq2 handles library size differences through its internal normalization procedure (McMurdie & Holmes, 2014). Prior to analysis, taxa were filtered to retain only species belonging to genera represented in the Unified Human Gastrointestinal Genome catalogue (Almeida et al., 2021), removing Kraken2 false positives arising from environmental organisms in the reference database. A design formula of `~ diet_group` was used with Omnivore as the reference level.

DESeq2 was run with Cook's distance outlier filtering disabled (`cooksCutoff = FALSE`) to recover taxa flagged as outliers due to the small sample size, which would otherwise remove potentially real biological signals. Log2 fold changes were subsequently shrunk using the apeglm method (Zhu et al., 2018) to stabilize estimates for low-count taxa and reduce inflation of fold changes from sparse counts. Wald test p-values were adjusted for multiple comparisons using the Benjamini-Hochberg procedure. Taxa with adjusted p-value < 0.05 and |log2 fold change| > 1 after shrinkage were considered differentially abundant. Results were visualized as an MA plot and a ranked bar plot of significant taxa. All code is provided in `scripts/05_deseq2.R`.

## References

- Almeida, A., et al. (2021). A unified catalog of 204,938 reference genomes from the human gut microbiome. *Nature Biotechnology*, 39, 105–114.
- Blanco-Míguez, A., et al. (2023). Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nature Biotechnology*, 41, 1633–1644.
- Bolger, A.M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics*, 30(15), 2114–2120.
- Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884–i890.
- De Filippis, F., Pasolli, E., Tett, A., et al. (2019). Distinct genetic and functional traits of human intestinal *Prevotella copri* strains are associated with different habitual diets. *Cell Host & Microbe*, 25(3), 444–453.
- Fackelmann, G., Manghi, P., Carlino, N., et al. (2025). Gut microbiome signatures of vegan, vegetarian and omnivore diets and associated health outcomes across 21,561 individuals. *Nature Microbiology*, 10, 40–55.
- Fernandes, A.D., et al. (2014). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome*, 2, 15.
- Lin, H., & Peddada, S.D. (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*, 11, 3514.
- Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*, 15, 550.
- Lu, J., Breitwieser, F.P., Thielen, P., & Salzberg, S.L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- Mallick, H., et al. (2021). Multivariable association discovery in population-scale meta-omics studies. *PLOS Computational Biology*, 17(11), e1009442.
- McMurdie, P.J., & Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE*, 8(4), e61217.
- McMurdie, P.J., & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. *PLOS Computational Biology*, 10(4), e1003531.
- Oksanen, J., et al. (2022). vegan: Community Ecology Package. R package version 2.6-2.
- Quince, C., Walker, A.W., Simpson, J.T., Loman, N.J., & Segata, N. (2017). Shotgun metagenomics, from sampling to analysis. *Nature Biotechnology*, 35, 833–844.
- Sonnenburg, J.L., & Bäckhed, F. (2016). Diet–microbiota interactions as moderators of human metabolism. *Nature*, 535, 56–64.
- Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York.
- Wood, D.E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
- Zhu, A., Ibrahim, J.G., & Love, M.I. (2018). Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics*, 35(12), 2084–2092.
