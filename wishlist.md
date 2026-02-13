# Pilot Result Wishlist (updated 13 Feb 2026)

*Prioritised by impact on the application × feasibility before 19 Feb deadline.*

## Already in the draft (confirmed results)
- ✅ MAGMA enrichments: neuron+spermatid p=0.030, XCI escape p=0.027, gametologs p=0.026
- ✅ Relate: 9/35 selected genes are ASD-associated, p=5.95×10⁻⁵; spermatid+brain vs brain-only p=0.002
- ✅ SFARI ASD genes enriched among X/Y spermatid DEGs, p=0.004
- ✅ X/Y spermatid DEGs enriched at A/B compartment borders, p=0.003
- ✅ ECH90 selective sweeps overlap A/B borders in spermatogenesis stages, p=0.005–0.028
- ✅ Y haplogroup stratification (I vs R) affects GWAS associations (qualitative)
- ✅ Triple-overlap Venn diagram: 5 genes (CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1) in all three sets (SFARI × selection × X/Y spermatid DEGs). All pairwise overlaps significant (***,***,**). Figure produced.

## Tier 1 — High impact, likely feasible before deadline

### 1. ✅ DONE — Triple-overlap Venn diagram (for Figure 2)
Result: 5 genes in triple overlap: CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1. All pairwise overlaps significant. Figure produced. Better than expected — 5 genes, all well-known ASD genes.

### 2. Expression trajectory plots for top candidates (for Figure 2)
**Data needed:** scRNA-seq expression values across spermatogenesis stages for candidate genes
**Analysis:** Plot expression of e.g. CASK, HUWE1, CDKL5, FRMPD4, DDX3X across spermatogonia → leptotene → pachytene → round spermatid → spermatozoa
**Why it matters:** If these genes show the signature pattern (expression → silencing at pachytene → reactivation in spermatids), that is visual proof of MSCI escape. Intuitive even for non-specialists.
**Effort:** Low-moderate — data is processed, need to extract and plot specific genes.

### 3. Haplogroup-specific hits near euchromatin–heterochromatin boundaries (strengthens WP2)
**Data needed:** Y-stratified GWAS results (already have), T2T-CHM13 pericentromeric boundary annotations
**Analysis:** Test whether haplogroup I-specific significant SNPs are enriched near pericentromeric boundaries compared to haplogroup R-specific SNPs.
**Why it matters:** Directly tests the sink model's key prediction. A significant enrichment would transform WP2 from "promising idea" to "supported hypothesis."
**Effort:** Moderate — need to define boundaries from T2T and run enrichment on existing GWAS output.

### 4. Haplogroup-specific hits in H3K9me3 brain regions (strengthens WP2)
**Data needed:** Y-stratified GWAS results, Roadmap Epigenomics H3K9me3 annotations for brain tissues
**Analysis:** Test whether haplogroup I-specific hits map to regions normally marked by H3K9me3 in brain.
**Why it matters:** The sink model predicts that haplogroup I (more heterochromatin → more sequestration of HP1/SUV39H) leads to "leaky" heterochromatin. SNPs in normally silent H3K9me3 regions reaching significance only in haplogroup I background would be strong evidence.
**Effort:** Moderate — Roadmap annotations are public, need to intersect with GWAS output.

## Tier 2 — High impact, moderate effort

### 5. X–Y interaction in iPSYCH
**Data needed:** X chromosome genotypes + Y haplogroup assignments for iPSYCH males
**Analysis:** Test whether effect sizes of X-linked ASD risk variants differ by Y haplogroup. Could use GEM or simple interaction model.
**Why it matters:** Would be the first evidence that both sex chromosomes jointly modulate ASD risk. Even a suggestive signal would justify WP3.3 and make the "unified framework" argument concrete.
**Effort:** Moderate-high — requires X chromosome GWAS to be run first (or use existing X-linked candidate SNPs).

### 6. X chromosome PRS adds to autosomal PRS
**Data needed:** iPSYCH X chromosome genotypes, autosomal PRS already computed
**Analysis:** Build X-linked PRS, test incremental R² beyond autosomal PRS.
**Why it matters:** Directly quantifies what the field is missing by excluding X from GWAS. Even a small but significant ΔR² makes a powerful argument for WP1.1.
**Effort:** Moderate — depends on how far Shannon's X chromosome GWAS has progressed.

### 7. Female XCI-escape enrichment
**Data needed:** iPSYCH female cases/controls, XCI escape gene annotations (from Tukiainen et al. 2017 or similar)
**Analysis:** In females, test whether ASD association is stronger among X-linked genes that escape XCI vs those that don't.
**Why it matters:** Directly tests the dosage-conflict prediction from a different angle — if selfish genes escape both MSCI and XCI, females get excess dosage from escaping genes. Already partially supported by the MAGMA XCI-escape result (p=0.027) but a female-specific test would be a clean replication.
**Effort:** Moderate — depends on female GWAS status.

### 8. Temporal dynamics of selection
**Data needed:** Existing Relate genealogies from intern thesis
**Analysis:** For the 9 ASD-associated positively selected genes, estimate timing of selection episodes. Do they cluster at the 45–55 kBP burst identified by Skov et al.?
**Why it matters:** Would show that the same sweep episode that displaced Neanderthal admixture also drove ASD-relevant genes to high frequency — a dramatic connection between archaic admixture and modern disease risk.
**Effort:** Moderate — Relate trees exist, need CLUES or similar for timing estimation.

## Tier 3 — Would be nice, harder before deadline

### 9. Autosomal chromatin differences between X- and Y-bearing spermatozoa
**Data needed:** Hi-C or ATAC-seq from FACS-sorted X- vs Y-bearing sperm (may not exist)
**Analysis:** Compare A/B compartment profiles on autosomes between the two cell types.
**Why it matters:** Direct evidence that Y chromosome presence/absence alters autosomal chromatin state — the most direct test of the sink model in humans.
**Effort:** High — depends on whether sorted sperm Hi-C data exists or can be generated.

### 10. Population differences in X-linked ASD architecture
**Data needed:** Population-stratified GWAS or 1000 Genomes data
**Analysis:** Test whether populations with different X chromosome selection histories have different X-linked ASD genetic architectures.
**Why it matters:** Supports the evolutionary dynamics argument in WP3.
**Effort:** High — requires multi-population analysis.

### 11. Doublecortin family deep dive
**Data needed:** Selection signals (Relate), expression data (scRNA-seq), ASD association (SFARI/iPSYCH) for DCX, DCLK1, DCLK2, DCX2
**Analysis:** Show that the DCX family sits at the intersection of neuronal migration, sperm motility, recent selection, and ASD association.
**Why it matters:** A compelling case study for the narrative and a potential figure panel. DCX regulates both microtubule dynamics in growing neurons and flagellar beating in sperm.
**Effort:** Moderate — mostly synthesis of existing results across the gene family.
