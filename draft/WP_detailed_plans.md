# Work Package Detailed Plans

Working document — not for submission. Purpose: nail down the concrete steps for each WP so the project description prose is grounded in reality.

**Personnel key:** - PD1 (Jan 2027 – Dec 2029): Postdoc 1 - PD2 (Jul 2027 – Jun 2030): Postdoc 2 - PD3 (Jan 2028 – Dec 2030): Postdoc 3 - PhD1 (Aug 2028 – Jul 2032): PhD student 1 - PhD2 (Aug 2029 – Jul 2032): PhD student 2

------------------------------------------------------------------------

Comparative genomics, Chromatin structure (HiC), GWAS, Access to protected cohorts, Primate genome access, single-cell expression

For the postdocs and phd students to reap the fruits of the projects they lead. They cannot run in sequence. They must run side-by-side with the projects they depend on. The postdoc expertise in each area must be present in the group at the same time.

1. GWAS
   1. Sex-stratified GWAS and XWAS in iPSYCH.
   2. Stratify by Y haplogroup identified from diagnostic SNPs in iPSYCH panel.
   3. Quantify missing heritability accounted for by epistatic effect of Y using a Bayesian model including all haplogroups and accounting for their genealogical relationship.
   4. Compute ASD genetic correlation between I and R haplogroup males.
   5. [Associate Y haplotype specific GWAS hits with genes at A/B compartments only exposed in either X or Y spermatids.]
   6. [Estimate the population-averaged polygenic risk score explained by hitchhiking-elevated frequencies of risk variants.]
   7. Repeat study on UK-biobank using fully sequenced Y chromosomes.
   8. [Develop method to predict polygenic risk score from autosomal and X variants on the background of known Y chromosome repeat profile.]
2. X chromosome drive
   1. Selection scans for driver genes using RELATE/SINGER
   2. Replicate scans in baboons
   3. Characterize dual roles in post-meiotic spermatogenesis and neurodevelopment for each driver candidate.
   4. Develop genealogical statistics sensitive to transient selection and robust to population admixture.
   5. [Estimate change in allele frequency by hitchhiking for each ASD associated variant.]
   6. [Compute GWAS effect size stratified by population frequency to estimate how much the hitchhiking elevates the population averaged polygenic risk score.]
3. Chromatin structure
   1. Reanalysis of single-cell Hi-C sperm data
   2. Bayesian inference of A/B borders
   3. Identification of affected TADs at A/B borders
   4. Identification of superloop structures spanning compartments on X
4. Y chromosome repeat content and prediction


## WP1.1: X Chromosome GWAS in Clinical Cohorts

**Lead:** PD1 (core work years 1–3) **Timeline:** 2027–2029

### What exactly needs to happen

1.  **Data access and QC (months 1–4)**
    -   Obtain iPSYCH genotype data for \~17,000 male and \~5,500 female ASD cases + matched controls
    -   X chromosome SNPs are typically excluded from standard QC pipelines — need custom QC: check for heterozygous calls in males (genotyping error), sex-discordant samples, pseudoautosomal region handling, HWE testing in females only
    -   Impute X chromosome using the TOPMed or 1000 Genomes Phase 3 reference panel (separate from autosomal imputation since X requires hemizygous calling in males)
    -   Decide dosage model: code male genotypes as 0/2 (full dosage compensation assumed) or 0/1 (no compensation), or test both and compare
2.  **Ancestry modelling on the X (months 3–8)**
    -   Standard GWAS controls for population structure using autosomal PCs, but the X has its own ancestry history (different Ne, different recombination, potential sex-biased migration)
    -   Implement the Zhang et al. 2023 approach: embed linear mixed model within a population genomic model that directly estimates X-linked ancestry
    -   Validate by comparing X-specific PCs vs autosomal PCs — if they diverge, ancestry correction matters
    -   Compute X-specific kinship matrix for the LMM
3.  **Sex-stratified GWAS (months 6–12)**
    -   Run GWAS separately in males and females
    -   Males: test each SNP under hemizygous model, adjusting for X-specific PCs and covariates (age, batch, array)
    -   Females: standard additive model but must account for XCI — compare dosage-compensated vs. escape-corrected models
    -   Meta-analyse male + female results with sex-specific effect sizes
    -   Generate Manhattan and QQ plots; identify genome-wide significant loci and suggestive loci
4.  **Heritability partitioning (months 8–14)**
    -   Run LDSC-SEG on X chromosome, partitioning heritability into:
        -   Spermatid/neuron overlap genes (from existing gene sets)
        -   XCI escapee genes (Tukiainen et al. 2017 annotations)
        -   Gametolog genes (X genes with Y homologs)
        -   ECH90 selective sweep regions (Anonymized Ref 1)
    -   Compare X-linked heritability enrichment per base pair against autosomal average
    -   Key question: is X-linked heritability enriched in spermatid-expressed gene sets specifically?
5.  **X-linked PRS (months 12–18)**
    -   Build X-linked polygenic risk score using PRSice-2 or PRS-CS adapted for X
    -   Test incremental R² of X-PRS beyond autosomal PRS in held-out iPSYCH sample
    -   Test sex-specific prediction: does X-PRS predict differently in males vs females?
6.  **Hitchhiking test (months 14–20)**
    -   For each X-linked SNP, compute selection coefficient using RELATE genealogies + CLUES
    -   Correlate selection strength with GWAS effect size and risk-allele frequency
    -   Positive correlation = hitchhiking prediction confirmed
    -   Control: repeat on autosomes (expect no correlation)
7.  **Multi-cohort pooling (months 18–30)**
    -   Obtain X chromosome summary statistics from PGC ASD GWAS and SPARK
    -   Meta-analyse with iPSYCH results using METAL
    -   Population diversity: PGC includes European, East Asian, and admixed samples

### Key outputs

-   First large-scale X-chr GWAS for ASD (publication)
-   X-linked heritability estimates partitioned by gene function
-   X-linked PRS with quantified incremental value
-   Selection–association correlation (hitchhiking evidence)

### Dependencies

-   Needs iPSYCH data access (already established)
-   PGC/SPARK data sharing agreements (standard, but takes time — initiate in month 1)

### Risks

-   X-chr imputation quality may be lower than autosomal (mitigate: use multiple reference panels, strict info score filtering)
-   Power: X has \~1,100 genes vs \~20,000 autosomal — individual loci may not reach significance, but gene-set tests should be well-powered with 17k cases

------------------------------------------------------------------------

## WP1.2: Identifying MSCI Escapees, Detecting Meiotic Drive, and Comparative Genomics

**Lead:** PD1 (years 1–3), with PD3 contributing comparative genomics in year 2 **Timeline:** 2027–2030

### What exactly needs to happen

1.  **Identify MSCI escapees from scRNA-seq (months 1–8)**
    -   Source data: published scRNA-seq from human and mouse spermatogenesis (e.g., Hermann et al. 2018, Shami et al. 2020)
    -   For each X-linked gene, extract expression across spermatogenesis stages: spermatogonia → leptotene → zygotene → pachytene → round spermatid → elongating spermatid → spermatozoa
    -   Classify genes as: (a) silenced at pachytene and remaining silent (MSCI maintained), (b) silenced at pachytene but reactivated in spermatids (MSCI escapees), (c) never silenced
    -   Category (b) are the meiotic drive candidates — post-meiotic expression is required for a gene to influence its own transmission
2.  **Expression trajectory plots for top candidates (months 4–8)**
    -   Plot trajectories for the 5 triple-overlap genes (CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1) plus other candidates (CASK, DDX3X, FRMPD4, DCX)
    -   Visual: x-axis = spermatogenesis stage, y-axis = normalised expression, shading for MSCI window
    -   If these show the silencing→reactivation pattern, that's visual proof of MSCI escape → Figure panel
3.  **Enrichment testing (months 6–12)**
    -   Test whether MSCI escapees are enriched for: SFARI ASD genes, iPSYCH GWAS hits (from WP1.1), positively selected genes (Relate), genes at A/B compartment borders
    -   Fisher's exact tests + permutation-based enrichment (shuffling gene labels while preserving chromosomal structure)
4.  **Functional characterization of driver candidates (months 8–24)**
    -   For the 5 triple-overlap genes, investigate known biology:
        -   **Cytoplasmic bridge transport:** Spermatids are connected by intercellular bridges (TEX14-dependent). Transcripts and proteins can pass between X- and Y-bearing spermatids. A driver gene would need to either (a) not pass through bridges (cell-autonomous) or (b) bias transport. Look for evidence of directed RNA transport, asymmetric protein localisation in published proteomics/imaging data
        -   **Mitochondrial asymmetry:** X- and Y-bearing spermatids may receive unequal mitochondria through biased bridge transport. If driver genes influence mitochondrial distribution (e.g. through cytoskeletal regulation), this could create transmission distortion via differential motility
        -   **Chromatoid body storage:** The chromatoid body is a phase-separated RNA granule in spermatids. Differentially stored mRNAs may be sequestered from translation in one spermatid type but not the other. Check whether candidate driver genes are found in chromatoid body transcriptomes (e.g., Meikar et al. 2014 datasets)
    -   This is largely bioinformatic synthesis of published data, not new experiments
5.  **Comparative selection scans in primates (months 12–30)**
    -   **Data:** publicly available great ape genomes (human, chimpanzee, gorilla, orangutan from the Great Ape Genome Project); baboon genomes from consortium participation (Baboon Genome Diversity Project or similar)
    -   **Method:** Run Relate on X chromosome alignments for each primate species separately, then compare positively selected gene sets across species
    -   **Key test:** Are the same X-linked gene sets (especially the 5 triple-overlap genes) independently under positive selection in macaques and/or baboons? Convergent selection across species that diverged 25–30 Mya would be strong evidence for intragenomic conflict (ecology changed, conflict persists)
    -   Control: compare X vs autosomal convergence rates — X should show excess convergence if conflict is the driver
    -   **Macaque half-brother analysis:** Reproduce/extend the finding that paternal genotype has an outsized contribution to heritability of sub-sociality traits. Map the specific paternal genomic regions contributing, focusing on sex chromosomes
6.  **Genealogy-based drive detection (months 18–30)**
    -   Apply ARG inference (Relate, tsinfer/tsdate, ARGweaver) to X chromosome data from 1000 Genomes + additional population panels
    -   Detect ongoing selective sweeps: regions where genealogies are shallower than expected (recent common ancestor = recent sweep)
    -   Compare drive signatures across populations (European, East Asian, African) — if drive is ongoing, sweep signatures should be population-specific; if ancient, they should be shared
    -   Cross-reference with ASD association: are ongoing sweeps enriched near ASD risk loci?

### Key outputs

-   Catalogue of MSCI escapees with expression trajectory evidence
-   Functional dossier on the 5 triple-overlap driver candidates
-   Cross-primate selection scan (publication)
-   Evidence for ongoing vs ancient drive on the human X

### Dependencies

-   Baboon genome access through consortium (confirm data sharing)
-   WP1.1 GWAS hits for enrichment testing

------------------------------------------------------------------------

## WP1.3: The Neuron/Spermatid Conflict — Autosomal Impact and Compensatory Evolution

**Lead:** PD3 (years 2–3), PhD1 contributes from year 3 **Timeline:** 2028–2030

### What exactly needs to happen

1.  **Autosomal selection scan on 1000 Genomes (months 1–6)**
    -   Run Relate on all autosomes using 1000 Genomes Phase 3 data (2,504 individuals, 26 populations)
    -   Identify genes under recent positive selection on autosomes
    -   Already done for X (intern thesis) — extend to autosomes using same pipeline and parameters for comparability
2.  **Compensatory evolution test (months 4–12)**
    -   Define "functionally linked to X/Y" gene sets:
        -   PPI partners of X-linked driver candidates (from STRING, BioGRID)
        -   Genes in the same pathways as driver candidates (KEGG, Reactome)
        -   Genes co-expressed with driver candidates in spermatogenesis (from scRNA-seq co-expression networks)
    -   Test: are autosomal PPI partners of X-linked driver genes enriched for positive selection compared to autosomal genes of matched size and expression level?
    -   This tests the prediction that autosomal genes continually compensate for transmission distortion caused by X-linked drivers — an arms race signature
3.  **Pathogenicity estimation (months 6–14)**
    -   For all genes in the spermatid/neuron overlap set, compute:
        -   CADD scores for segregating variants
        -   pLI (probability of loss-of-function intolerance) scores
        -   Missense constraint (from gnomAD)
    -   Test: do spermatogenesis-involved genes accumulate more pathogenic variants than expected? The prediction: genes under meiotic drive have inflated frequencies of linked deleterious variants (hitchhiking)
    -   Compare pathogenic variant burden in the spermatid/neuron overlap set vs. neuron-only or spermatid-only genes
4.  **In-frame stop codon analysis (months 10–16)**
    -   For autosomal genes identified as potential co-optees of X-linked drivers:
        -   Count in-frame stop codons from gnomAD
        -   Test whether these genes lose function more often (convergent pseudogenisation, like TTLL11)
    -   Rationale: if autosomal genes are co-opted by drivers and this co-option is deleterious, selection may favour loss-of-function — analogous to TTLL11 inactivation
5.  **Doublecortin (DCX) family case study (months 12–20)**
    -   Compile evidence across DCX, DCLK1, DCLK2:
        -   Selection signals (from Relate scan)
        -   Expression across spermatogenesis stages and brain regions
        -   ASD association (SFARI, iPSYCH GWAS)
        -   Protein function: microtubule stabilisation in neurons + flagellar beating in sperm
    -   DCX is X-linked; DCLK1 and DCLK2 are autosomal paralogs — test whether the autosomal paralogs show compensatory selection patterns
6.  **TTLL11 as a defence model (months 14–22)**
    -   Synthesise evidence for TTLL11 convergent inactivation in humans + gorillas
    -   Map TTLL11 interactors and substrates — which tubulin isoforms does it modify?
    -   Test whether TTLL11 substrates overlap with bridge-transport genes or driver candidate interactors
    -   This provides a concrete mechanistic narrative: drive gene co-opts tubulin transport → host inactivates the tubulin modifier → defense
7.  **Multivariate PGS and ASD subtyping (months 18–28)**
    -   Build PGS for the meiotic drive gene set (X-linked spermatid/neuron overlap genes)
    -   Test whether this PGS explains phenotypic variation across ASD subtypes (e.g., language delay, intellectual disability, social communication difficulty)
    -   Use iPSYCH phenotype data — ASD subtypes defined by ICD-10 codes
    -   Regression: ASD severity/subtype \~ autosomal PGS + drive-gene PGS + covariates

### Key outputs

-   Autosomal selection scan + compensatory evolution test (publication)
-   Pathogenicity burden comparison across gene sets
-   DCX family case study (narrative for project description)
-   Drive-gene PGS and ASD subtype association

### Dependencies

-   WP1.1 GWAS results for PGS construction
-   WP1.2 candidate gene lists for enrichment base

------------------------------------------------------------------------

## WP1.4: 3D Chromatin Architecture and Meiotic Drive

**Lead:** PD3 (primary, years 2–3), PhD1 (Part A, years 2–3) **Timeline:** 2028–2030

### What exactly needs to happen

1.  **Obtain and process Hi-C data (months 1–4)**
    -   Source: published Hi-C from macaque spermatogenesis stages (Alavattam et al. 2019 or similar), human sperm Hi-C (Ke et al. 2017), possibly human spermatogenesis Hi-C as new datasets emerge
    -   Process with HiC-Pro or Juicer → generate contact matrices at multiple resolutions (10kb, 50kb, 100kb)
    -   Call A/B compartments using eigendecomposition of the correlation matrix (standard Lieberman-Aiden method)
2.  **Bayesian compartment border detection (months 3–10)**
    -   Apply the Bayesian border detection method developed in Anonymized Reference 5 (master thesis)
    -   This method provides posterior probabilities for border positions rather than sharp calls, enabling fine-scale mapping
    -   Run on each spermatogenesis stage separately: spermatogonia, leptotene/zygotene, pachytene spermatocytes, round spermatids, elongating spermatids, spermatozoa
    -   Generate border probability profiles along the X chromosome for each stage
3.  **Characterize fluid vs stable borders (months 6–14)**
    -   Key analysis from KM's v3 comments: which A/B borders are shared across tissues (constitutive) and which change during spermatogenesis (dynamic)?
    -   Compare spermatogenesis-stage borders to borders in:
        -   Somatic tissues (brain, liver, lung — from Roadmap Epigenomics Hi-C)
        -   Female cells (to distinguish sex-specific from spermatogenesis-specific)
    -   Classify borders as: constitutive (same across all tissues/stages), spermatogenesis-specific (present only during meiosis), or stage-specific (present only in particular spermatogenesis stages)
    -   Prediction: driver genes should cluster at spermatogenesis-specific borders — these are the borders that create the opportunity for MSCI escape
4.  **Test enrichment of driver candidates at borders (months 10–18)**
    -   For each border class (constitutive, spermatogenesis-specific, stage-specific):
        -   Test overlap with MSCI escapees (from WP1.2)
        -   Test overlap with positively selected genes (from WP1.2 Relate scan)
        -   Test overlap with SFARI ASD genes
        -   Test overlap with the 5 triple-overlap driver candidates
    -   Permutation testing: shuffle gene positions along the X while preserving gene density to compute empirical p-values
    -   The key prediction: MSCI escapees should be enriched at spermatogenesis-specific borders, not constitutive ones
5.  **A/B eigenvector figure (months 12–18)**
    -   Produce the figure from KM's Potential Figures list: A/B eigenvector plot showing X chromosome compartments in sperm, aligned with chromosome ideogram, overlaid with positions of positively selected genes
    -   Show how selection targets coincide with compartment boundaries
    -   This is both a result and a key visual for the application/publications
6.  **Extend to human spermatogenesis data (months 18–30)**
    -   As new human spermatogenesis Hi-C datasets become available (several groups are generating these), repeat the analysis
    -   Compare human and macaque border profiles — conserved borders suggest functional importance
    -   Compare compartment dynamics between X and autosomes during meiosis — X should show more dramatic reorganisation due to MSCI

### Key outputs

-   Classification of A/B borders by stability across tissues/stages
-   Enrichment of driver candidates at dynamic borders (key mechanistic result)
-   A/B eigenvector figure for publications
-   Bayesian border detection method (software release)

### Dependencies

-   Builds on existing macaque Hi-C data + border method (Anonymized Ref 5)
-   Benefits from WP1.2 gene lists

------------------------------------------------------------------------

## WP2.1: Y Haplogroup-Stratified GWAS and Chromatin-State Mapping

**Lead:** PD2 (Jul 2027 – Jun 2030) **Timeline:** 2027–2030 (ramps up from Jul 2027)

### What exactly needs to happen

1.  **Y haplogroup assignment (months 1–4, from Jul 2027)**
    -   Extract Y-chromosomal SNPs from iPSYCH genotyping data (males only, \~17,000 ASD + controls)
    -   Assign haplogroup using yhaplo, Y-LineageTracker, or similar tool
    -   Validate against the T2T-CHM13 Y chromosome assemblies — 43 complete Y assemblies from Rhie et al. 2023 provide ground-truth haplogroup assignments with full structural characterisation
    -   Classify into I, R (and subgroups I1, I2, R1a, R1b) — these two dominate in Danish males
    -   For rare haplogroups (E, J, G, N): group as "other" or exclude from primary analysis
2.  **Haplogroup-stratified autosomal GWAS (months 4–12)**
    -   Split the male cohort by haplogroup (I vs R as primary comparison)
    -   Run standard autosomal GWAS within each stratum using PLINK2
    -   Key comparison: which SNPs are significant in haplogroup I but not R, and vice versa?
    -   Compute heterogeneity statistics (Cochran's Q) to identify SNPs with significantly different effects between haplogroups
3.  **Interaction testing with GEM (months 8–14)**
    -   Run genome-wide SNP × haplogroup interaction using GEM (Gene-Environment interaction in Millions, Wang et al. 2021)
    -   Treat haplogroup as the "environment" — models the interaction between each autosomal SNP and Y haplogroup on ASD risk
    -   GEM is designed for biobank-scale data and handles related individuals efficiently
    -   Identify SNPs with significant interaction terms (FDR \< 0.05 genome-wide)
4.  **A/B compartment mapping by Y haplotype (months 10–20)**
    -   This is the direct mechanistic test of the sink model (from KM's v3 comments)
    -   Approach: use expression data (GTEx or iPSYCH if eQTL data available) to infer chromatin compartment status indirectly
    -   More direct approach: if ATAC-seq or Hi-C data from males with known Y haplogroup becomes available, compare A/B compartments between I and R males
    -   Alternatively: use the haplogroup-stratified GWAS results — regions where haplogroup I shows association but R does not may correspond to regions where chromatin is "opened" by the I haplogroup's greater heterochromatin sink effect
    -   Map these regions to T2T-CHM13 annotations to check if they lie near pericentromeric boundaries
5.  **Bayesian heritability estimation stratified by Y haplotype (months 14–22)**
    -   Estimate SNP-heritability of ASD separately for haplogroup I and R males
    -   Use GCTA-GREML or LDSC, computing GRM/LD scores within each stratum
    -   Key question: is ASD more heritable in one haplogroup background than the other? The sink model predicts that haplogroup I males have more "exposed" risk variants → potentially higher heritability
    -   Also partition heritability by functional annotation within each haplogroup
6.  **Extension to other cohorts (months 20–30)**
    -   UK Biobank: assign Y haplogroups (already available in the UK Biobank Y-chr analysis), test haplogroup × autosomal SNP interactions for ASD-related traits (social communication, repetitive behaviours from questionnaire data)
    -   African-ancestry cohorts (if available through PGC): deeply diverged Y haplogroups (A, B, E) provide a more extreme test — much greater heterochromatin divergence between haplogroups

### Key outputs

-   Y-stratified GWAS of ASD (publication — first of its kind)
-   SNP × haplogroup interaction results
-   Heritability estimates stratified by Y haplogroup
-   Identification of autosomal regions with haplogroup-dependent effects

### Dependencies

-   iPSYCH data access (shared with WP1.1 — PD1 and PD2 coordinate data prep)
-   T2T-CHM13 annotations (public)

------------------------------------------------------------------------

## WP2.2: Chromatin Landscape Analysis

**Lead:** PD2 (years 2–3) **Timeline:** 2028–2030

### What exactly needs to happen

1.  **Map haplogroup-specific SNPs to chromatin states (months 1–6)**
    -   Download Roadmap Epigenomics 15-state ChromHMM annotations for brain tissues: prefrontal cortex (E073), hippocampus (E071), cerebellum (E068), anterior caudate (E070)
    -   For each haplogroup-specific significant SNP (from WP2.1), annotate with chromatin state in each brain tissue
    -   Compare chromatin state distribution of haplogroup I-specific vs R-specific SNPs
2.  **H3K9me3 enrichment test (months 4–10)**
    -   The key sink model prediction: haplogroup I (more heterochromatin → more HP1/SUV39H sequestration) leads to "leaky" heterochromatin
    -   Test: are haplogroup I-specific SNPs enriched in regions normally marked by H3K9me3 in brain?
    -   Use LDSC-SEG with H3K9me3 ChIP-seq annotations from Roadmap/ENCODE
    -   Also test GoShifter for locus-level enrichment
    -   Compare against H3K4me3 (active marks) as control — the sink model predicts effects specific to heterochromatin, not euchromatin
3.  **Pericentromeric boundary enrichment (months 6–14)**
    -   Use T2T-CHM13 with its complete centromeric and pericentromeric annotations
    -   Define euchromatin–heterochromatin boundaries from the T2T assembly: positions where α-satellite/satellite repeat density transitions to gene-containing euchromatin
    -   Test whether haplogroup I-specific SNPs are enriched near these boundaries (within 1 Mb, 500 kb, 100 kb — varying distance thresholds)
    -   Haplogroup R-specific SNPs should show no such enrichment (or weaker)
    -   This is the most direct genomic test of the sink model in humans
4.  **Differential methylation analysis (months 10–20, contingent on data)**
    -   IF methylation array data with Y haplogroup information is available (possible through iPSYCH or external collaborators):
        -   Assign Y haplogroups to methylation-typed males
        -   Run differential methylation analysis (DMRcate, minfi) comparing haplogroup I vs R
        -   Focus on: pericentromeric regions, ASD gene promoters, enhancers active in brain
    -   IF not available: this sub-aim becomes a future direction
5.  **Integration with WP1.4 (months 14–20)**
    -   Compare autosomal A/B compartment borders (from WP1.4 methods) with haplogroup-specific association signals
    -   Test: do haplogroup-specific SNPs concentrate at A/B compartment boundaries in brain tissue?
    -   This bridges the X-chromosome chromatin story (WP1.4) with the Y-chromosome epistasis story (WP2)

### Key outputs

-   Chromatin state annotation of haplogroup-specific associations
-   H3K9me3 enrichment test result (key for/against sink model)
-   Pericentromeric boundary enrichment (direct sink model test)

### Dependencies

-   WP2.1 haplogroup-stratified GWAS results (core dependency)
-   WP1.4 chromatin border methods (synergy)

------------------------------------------------------------------------

## WP2.3: Haplogroup-Specific Pathway and Network Analysis

**Lead:** PhD1 Part A (Aug 2028 – Jul 2030), PD2 contributes **Timeline:** 2029–2031

### What exactly needs to happen

1.  **Pathway mapping (months 1–6)**
    -   Take haplogroup-specific gene lists from WP2.1 (genes near significant SNPs in I-only, R-only, and shared categories)
    -   Map to GO terms, KEGG pathways, Reactome pathways using g:Profiler or clusterProfiler
    -   Focus areas: chromatin modification, synaptic function, neurodevelopmental pathways, immune regulation
2.  **Convergent pathway test (months 4–12)**
    -   Key question: do genes significant in haplogroup I and genes significant in haplogroup R converge on the same pathways despite being different genes?
    -   This would mean the same biology is being perturbed through different chromatin routes — strong evidence for the sink model
    -   Method: compute pathway overlap using hypergeometric test, then assess significance with permutation (shuffle gene-to-haplogroup assignments)
    -   Also test PPI network proximity: are I-specific and R-specific genes closer in PPI space than random gene sets of the same size?
3.  **Haplogroup-specific PRS (months 8–18)**
    -   Build separate PRS using haplogroup I-specific and R-specific effect sizes
    -   In held-out samples, compare prediction accuracy of:
        -   Standard (unstratified) autosomal PRS
        -   Haplogroup I-specific PRS (applied to I males only)
        -   Haplogroup R-specific PRS (applied to R males only)
    -   Key metric: incremental R² of stratified PRS over unstratified PRS
    -   This directly addresses "missing heritability": if stratified PRS outperforms unstratified, it means ignoring Y haplogroup loses heritable signal
4.  **Missing heritability quantification (months 14–22)**
    -   Estimate heritability using all males combined → H²_combined
    -   Estimate within haplogroup I → H²_I, within haplogroup R → H²_R
    -   Compute recovered heritability: (H²_I × N_I + H²_R × N_R) / N_total vs H²_combined
    -   If the sink model is correct, H²_stratified \> H²_combined because epistatic effects that cancel out in mixed analysis become additive within strata

### Key outputs

-   Pathway convergence analysis (publication)
-   Haplogroup-specific PRS with improved prediction
-   Quantification of recovered heritability

### Dependencies

-   Requires WP2.1 results as primary input
-   Benefits from WP2.2 chromatin annotations for interpretation

------------------------------------------------------------------------

## WP3.1: Population Genetic Models of the Neuron/Spermatid Conflict

**Lead:** PhD1 Part B (Aug 2030 – Jul 2032), PD3 contributes foundation in year 3 **Timeline:** 2029–2032

### What exactly needs to happen

1.  **Analytical model of X-linked meiotic drive with pleiotropy (months 1–8)**
    -   Develop a diploid population genetic model with:
        -   X-linked locus with alleles A (neutral) and D (driver)
        -   Transmission ratio distortion: D-bearing sperm have frequency k \> 0.5 (parameterise k from 0.5 to 1.0)
        -   Fitness cost: male carriers of D have reduced fitness w_m = 1 - s (neurodevelopmental cost), female carriers w_f depends on XCI escape
        -   Dosage-conflict: if D escapes XCI in females, females have excess expression → fitness cost in females too, but this drives compensatory reduction → male deficiency
    -   Derive equilibrium frequency of D as a function of k and s
    -   Compute steady-state ASD risk: given D frequency and penetrance, what fraction of males in the population carry ASD-relevant alleles maintained by drive?
2.  **Multi-locus extension (months 6–14)**
    -   Extend from single-locus to multi-locus model (many driver genes simultaneously)
    -   Parameterise using observed data: \~35 X-linked genes under selection, 9 of which are ASD-associated
    -   Each locus has drive strength k_i and fitness cost s_i
    -   Predict total polygenic ASD risk from the sum of driver-maintained alleles
    -   Compare predicted polygenic architecture with observed GWAS results from WP1.1
3.  **Hitchhiking model (months 10–18)**
    -   Model linked selection on the X: given a driver sweep, how much do nearby deleterious variants hitchhike?
    -   Compute the expected elevation in risk-allele frequency as a function of recombination distance from the driver
    -   Key empirical test (using WP1.1 data): estimate how many causal variants on X have hitchhiked to frequencies higher than their detrimental effect would predict
    -   This requires comparing observed allele frequencies with frequencies predicted under drift alone (using demographic models)
4.  **Quantify PGS from Y-dependent genes (months 14–22)**
    -   Using WP2.1 results: identify autosomal genes whose A/B status depends on Y haplotype
    -   Build a PGS restricted to these Y-dependent genes
    -   Test: does this Y-dependent PGS explain the difference in ASD prevalence between haplogroup I and R males?
    -   If yes: the sink model not only creates haplogroup-specific associations but quantitatively accounts for prevalence differences
5.  **Historical dynamics: the 45–55 kBP sweep (months 18–26)**
    -   Parameterise the model with the sweep episode from Anonymized Ref 1 (14 loci, 10% of X affected)
    -   Compute the expected increase in population-averaged polygenic ASD risk when these sweeps occurred
    -   Did the sweep episode that coincided with Neanderthal admixture introgression also inflate ASD-relevant variation?
    -   Temporal dynamics: plot ASD risk trajectory over last 100,000 years
6.  **Software release (months 24–30)**
    -   Package analytical models and simulation code as open-source tools
    -   Include parameter estimation from GWAS and selection data

### Key outputs

-   Formal theory of meiotic drive → ASD risk (publication — theoretical)
-   Quantitative predictions testable against data
-   Historical reconstruction of ASD risk dynamics
-   Open-source simulation framework

### Dependencies

-   WP1.1 GWAS results for parameterisation
-   WP1.2 selection estimates for drive strength
-   WP2.1 Y-haplotype-dependent gene sets

------------------------------------------------------------------------

## WP3.2: Y-Haplotype-Based Disease Variant Prediction

**Lead:** PhD2 (Aug 2029 – Jul 2032) **Timeline:** 2029–2032

### What exactly needs to happen

1.  **Build Y-haplotype-to-variant-exposure map (months 1–12)**
    -   Using WP2.1 results: for each Y haplogroup (I, R, and subgroups), compile the set of autosomal variants that show significant effects specifically in that haplogroup background
    -   Using WP2.2 results: identify the chromatin regions "exposed" (shifted from heterochromatic to euchromatic) in each haplogroup
    -   For each male, his Y haplotype predicts which risk variants are functionally active
    -   Build a predictive model: Y haplotype → predicted set of exposed variants → predicted ASD risk profile
2.  **Y repeat composition (months 6–16)**
    -   Beyond haplogroup assignment, incorporate Y repeat composition (DYZ1 length, α-satellite array size) using available structural Y chromosome data
    -   Test whether continuous Y repeat length is a better predictor than categorical haplogroup assignment
    -   This requires Y repeat estimation from array/sequencing data — methods exist (e.g., k-mer based estimation from short reads)
3.  **ASD subtype stratification (months 12–24)**
    -   Using the predicted variant exposure profiles, cluster ASD cases into subtypes
    -   Methods: k-means, hierarchical clustering, or latent class analysis on the predicted active-variant vectors
    -   Validate clusters against phenotypic data: do the Y-haplotype-defined subtypes differ in:
        -   Age of diagnosis
        -   Intellectual disability co-occurrence
        -   Language delay
        -   Severity measures
        -   Co-occurring conditions (ADHD, epilepsy, anxiety)
    -   If subtypes with distinct aetiologies emerge, this is the clinically translatable result
4.  **Pharmacogenomic implications (months 18–30)**
    -   For the identified subtypes, check whether the defining genetic pathways have known pharmacological targets:
        -   Glutamate signalling pathway → Lamotrigine
        -   Oxytocin signalling pathway → Oxytocin
        -   Chromatin modification pathway → HDAC inhibitors (preclinical)
    -   This is speculative but important for the "Future Applications" narrative
    -   Determine whether any clinical trials for ASD drugs show differential response that could be explained by Y haplogroup
5.  **Cross-validation and replication (months 24–34)**
    -   Train prediction model on iPSYCH, validate on UK Biobank (ASD-related traits) and/or SPARK
    -   Assess prediction accuracy: AUC for ASD case/control, correlation with symptom severity
    -   Compare Y-haplotype-informed prediction vs standard PRS

### Key outputs

-   Y-haplotype-to-variant-exposure mapping (primary deliverable)
-   ASD subtype stratification by Y-predicted variant profiles
-   Cross-validated prediction model
-   Discussion of pharmacogenomic implications

### Dependencies

-   Requires WP2.1 and WP2.2 results as core inputs
-   Benefits from WP3.1 models for interpreting allele frequency dynamics

------------------------------------------------------------------------

## Cross-WP Synergies Summary

| Interaction | Mechanism | Timing |
|------------------------|------------------------|------------------------|
| PD1 ↔ PD2 | Shared iPSYCH data infrastructure, QC pipelines, GWAS methods | Years 1–3 (continuous) |
| PD1 → PD3 | X-chr selection results seed autosomal compensatory test | Year 2 (PD1 results → PD3 starts) |
| PD2 → PD3 | Y chromatin landscape methods inform X chromatin analysis | Years 2–3 |
| PD3 → PhD1 | PD3 establishes pathway analysis framework, PhD1 extends it | Year 3 (handover) |
| PD1 + PD2 → PhD1 | Both GWAS results feed PhD1's pathway convergence analysis | Year 3 |
| PD2 → PhD2 | Y-stratified results are PhD2's primary input | Years 3–4 (handover) |
| PhD1 → PhD2 | Modelling framework (WP3.1) informs variant prediction model | Years 4–5 |
| WP1.4 ↔ WP2.2 | Same Hi-C/chromatin methods applied to X and autosomal questions | Years 2–3 |
| WP1.1 → WP3.1 | GWAS effect sizes parameterise population genetic models | Years 3–4 |
| WP2.1 → WP3.2 | Haplogroup-specific variants are the inputs to prediction | Years 3–4 |