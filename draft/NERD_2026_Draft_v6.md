# NERD 2026 Application — Working Draft v6

**Status:** Updated to match restructured WP plan (WP_detailed_plans.md, "New plan" section). Written for Committee for Natural and Technical Sciences audience.
**Deadline:** 19 February 2026, 14:00 CET

## Timeline

| Work package / Year                                      | 2027 | 2028 | 2029 | 2030 | 2031 | 2032 |
|:---------------------------------------------------------|:----:|:----:|:----:|:----:|:----:|:----:|
| WP1 GWAS of X chromosome and Y-stratified autosomes      |  X   |  X   |  X   |  X   |  -   |  -   |
| WP2 Detection and characterisation of X-linked drive      |  X   |  X   |  X   |  X   |  -   |  -   |
| WP3 3D chromatin architecture in spermatogenesis          |  -   |  X   |  X   |  X   |  -   |  -   |
| WP4 Y chromosome repeat content and prediction            |  -   |  X   |  X   |  X   |  X   |  X   |
|                                                           |      |      |      |      |      |      |
| **Personnel**                                             |      |      |      |      |      |      |
| Postdoc 1 (WP1, WP2)                                     |  X   |  X   |  X   |  -   |  -   |  -   |
| Postdoc 2 (WP1, WP3)                                     |  ½   |  X   |  X   |  ½   |  -   |  -   |
| Postdoc 3 (WP2, WP3)                                     |  -   |  X   |  X   |  X   |  -   |  -   |
| PhD 1 (WP1, WP4) — 4 yr                                  |  -   |  ½   |  X   |  X   |  X   |  ½   |
| PhD 2 (WP4) — 4 yr                                       |  -   |  -   |  ½   |  X   |  X   |  ½   |

------------------------------------------------------------------------

## Potential Figures

1. Three way Venn diagram:
   - Genes subject to recent selection
   - Genes differentially expressed in X an Y spermatids
   - Genes linked to ASD
2. A/B eigenvector plot showing X chromosome A/B compartments in sperm aligned with chromosome ideogram overlaid with gene positions for genes with recent positive selection
3. Expression trajectory plots for top candidates
4. Project timeline (Gantt chart)



## 1. PROJECT TITLE (max 150 characters)

**Option A (137 chars):** Selfish sex chromosomes and the genetic basis of autism: how intragenomic conflict shapes human neurodevelopment

**Option B (148 chars):** Sex chromosome conflict and autism: how meiotic drive on X and chromatin regulation by Y shape the genetic basis of neurodevelopment

**Option C (126 chars):** The evolutionary genetics of autism: selfish sex chromosomes and their consequences for human neurodevelopment

------------------------------------------------------------------------

## 2. BRIEF PROJECT DESCRIPTION (max 2,000 characters)

Autism spectrum disorder (ASD) is four times more common in boys than in girls — one of the most robust findings in psychiatry — yet no satisfactory genetic explanation exists for this sex bias. This project proposes that both sex chromosomes actively drive ASD risk through distinct evolutionary mechanisms, and that recognizing these mechanisms will transform our understanding of how ASD-relevant variation is maintained in human populations.

The core idea draws on intragenomic conflict theory. In males, sperm cells carrying the X chromosome compete against those carrying the Y. An X-linked gene that is expressed after the cell division separating X- and Y-bearing sperm can bias its own transmission — a phenomenon called meiotic drive. Because roughly one in four neuron genes is also expressed in sperm cells, this transmission advantage comes at a cost: variants favoured during sperm competition may be mildly harmful for brain development. We present pilot evidence that genes at this intersection — simultaneously under positive natural selection, differentially expressed between X- and Y-bearing sperm, and associated with ASD — overlap far more than expected by chance.

On the Y chromosome, a separate mechanism may operate. Human Y chromosomes differ dramatically in the amount of repetitive heterochromatin they carry. This structural variation may titrate chromatin-modifying proteins away from the rest of the genome, altering which genes are accessible and thereby creating Y-lineage-dependent differences in ASD risk.

The project combines population genomics, clinical cohort analysis of more than 22,000 ASD cases, comparative primate genomics, and formal evolutionary modelling.

------------------------------------------------------------------------

## 3. PROJECT DESCRIPTION (max 30,000 characters)

<!-- BEGIN PROJECT DESCRIPTION — prose below this line counts toward the 30,000-character limit -->

### The Problem

Autism spectrum disorder (ASD) affects nearly 3% of children in Denmark and shows a 4:1 ratio of affected boys to girls. This sex bias is among the most replicated observations in human genetics, yet genome-wide association studies (GWAS) — which scan the genome for variants correlated with disease risk — have largely excluded the sex chromosomes from analysis due to their non-standard inheritance. The X chromosome is present in two copies in females but only one in males; the Y chromosome is present only in males. As a result, standard statistical methods designed for autosomes (the 22 non-sex chromosomes) do not apply. This project develops the analytical framework needed to bring the sex chromosomes into the study of ASD, and tests two specific hypotheses about how evolutionary forces on X and Y contribute to autism risk.

### Hypothesis 1: Meiotic Drive on the X Chromosome

**The basic mechanism.** During the production of sperm cells (spermatogenesis), a precursor cell containing both X and Y divides to produce two types of haploid cells: sperm carrying the X, and sperm carrying the Y. Early in this process, a phenomenon called meiotic sex chromosome inactivation (MSCI) shuts down transcription from the sex chromosomes. MSCI is thought to prevent genes on the X from acting against Y-bearing cells or vice versa. However, some genes escape this silencing and are re-expressed in the haploid sperm cells (spermatids). A gene that escapes silencing in X-bearing spermatids can, in principle, bias its own transmission — for example by enhancing motility or survival of X-bearing sperm at the expense of Y-bearing sperm. This is meiotic drive: a violation of Mendel's law of equal segregation, driven by molecular competition within the germline.

Drive is widespread in nature. It has been documented in fungi, insects, and mammals, and theoretical models predict it should be especially potent on the X chromosome, which spends two-thirds of evolutionary time in females (where it recombines) and one-third in males (where it does not). The effective population size of the X is three-quarters that of autosomes, amplifying the effects of linked selection.

**Why drive is hard to observe directly.** A natural objection is that transmission distortion and sex-ratio bias have not been reported in humans. This is expected under the dynamics of drive. Each driver produces strong distortion only against particular Y chromosomes and on particular autosomal backgrounds. As the driver increases in frequency, the population's allele frequencies at interacting loci adjust, suppressing the average distortion. Drive is therefore transient: it flares, pushes the driver (and its hitchhiking neighbourhood) to high frequency, and then subsides as the genomic background equilibrates — leaving no ongoing sex-ratio bias but a lasting footprint of elevated allele frequencies and their collateral effects. This transient nature means drive is best detected through its evolutionary signatures (selective sweeps, hitchhiking) rather than through contemporary transmission ratios.

**The link to autism.** Proteomic surveys show that roughly 80% of brain proteins are also present in testis, and approximately one in four genes active in neurons is also expressed in spermatids [Matos et al. 2021]. This overlap creates an evolutionary conflict: a gene variant selected for its advantage in spermatogenesis may be mildly deleterious for brain development. Because the X chromosome carries a disproportionate fraction of neurodevelopmental genes — it accounts for approximately 20% of neuro-anatomical variation relevant to ASD despite constituting only 5% of the genome [Mallard et al. 2021] — this conflict falls disproportionately on the X.

Two specific predictions follow. First, the "hitchhiking" prediction: ASD risk variants physically close to driven genes on the X should reach higher population frequencies than variants with comparable effects on autosomes, because they are dragged along by the selective sweep of the driven allele. Second, the "dosage-conflict" prediction: because MSCI and somatic X chromosome inactivation (XCI, the process by which one X is silenced in female cells) share a conserved silencing mechanism [Hornecker et al. 2007], genes that escape MSCI in sperm will also tend to escape XCI in females. Since the X chromosome spends two-thirds of evolutionary time in females, selection will resolve the resulting female dosage excess by reducing expression. This leaves males — who have only one copy — with too little expression. For neurodevelopmental genes, this dosage deficiency creates a male-specific risk.

**Pilot evidence — five converging lines.** (i) MAGMA analysis of a large Danish clinical cohort shows that neuron genes co-expressed in spermatids are significantly enriched for ASD association in males (p = 0.030), while neuron genes without spermatid expression show no enrichment — confirming that specifically the dual-function genes carry risk [Anonymized Reference 2]. (ii) ASD association is elevated among genes that escape female XCI (p = 0.027), as the dosage-conflict model predicts [Anonymized Reference 2]. (iii) An analysis of evolutionary trees (genealogies) reconstructed from population genomic data identifies 35 X-linked genes under recent positive selection; 9 of these are ASD-associated, far exceeding chance (p = 5.95 × 10⁻⁵). Selection is significantly more frequent among genes expressed in both spermatids and brain than among brain-only genes (p = 0.002) [Anonymized Reference 3]. (iv) Known ASD genes are significantly enriched among genes differentially expressed between X-bearing and Y-bearing spermatids (p = 0.004) [Anonymized Reference 6]. Because post-meiotic spermatids share cytoplasm through intercellular bridges, most transcripts are equally distributed; genes that nonetheless differ must be regulated cell-autonomously — precisely the property needed for meiotic drive. (v) Intersecting three independent gene sets — ASD genes (SFARI database), genes under positive selection, and genes differentially expressed between X- and Y-bearing spermatids — yields five genes in the triple overlap (CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1), all of them well-established ASD risk genes [Anonymized References 3, 6]. Each pairwise and triple overlap is statistically significant.

**Convergent evidence from primate evolution.** The tubulin-modifying enzyme TTLL11 has been independently inactivated in both humans and gorillas, but remains functional in chimpanzees. TTLL11 polyglutamylates tubulin and is required for directed sperm motility; knockout mice produce sperm that swim in circles. One hypothesis is that TTLL11 inactivation represents a host defence: if an X-linked driver co-opted TTLL11-dependent microtubule transport to bias its own transmission, then organisms that destroyed the transport machinery would be freed from drive. The convergent loss in two independent great ape lineages — both with low sperm competition — is consistent with this scenario. Evidence from a macaque half-brother study further demonstrates that paternal genotype makes a disproportionate contribution to heritability of a sub-sociality trait — a phenotype correlated with ASD in human children — pointing directly at spermatogenesis as the arena in which heritable neurodevelopmental risk is shaped.

### Hypothesis 2: The Y Chromosome as a Chromatin Regulator

**The heterochromatin sink.** The second hypothesis concerns the Y chromosome not as a carrier of protein-coding genes — it has very few — but as a structural element that sequesters chromatin-modifying proteins. Large portions of the Y consist of repetitive DNA packaged as heterochromatin, a condensed form of chromatin associated with gene silencing. Maintaining heterochromatin requires specific proteins (notably HP1 and the methyltransferase SUV39H) that exist in limited quantities in the nucleus. A Y chromosome carrying more heterochromatin will sequester more of these proteins, leaving fewer available to maintain heterochromatin elsewhere in the genome. This "sink" effect has been rigorously demonstrated in Drosophila, where adding or removing Y chromosomes causes measurable genome-wide redistribution of histone marks at the boundaries between heterochromatin and euchromatin (the active form of chromatin) [Brown et al. 2020].

**Variation between human Y haplogroups.** Recent telomere-to-telomere assemblies of 43 human Y chromosomes [Rhie et al. 2023] reveal striking structural variation between Y lineages (haplogroups). Centromeric repeat arrays vary more than 2-fold in size between major haplogroups (e.g. haplogroup R1b: mean 341 kbp; haplogroup I: mean 787 kbp). The DYZ1 satellite repeat varies over an order of magnitude (7–98 Mb). In Denmark, most men carry a Y from either haplogroup I or haplogroup R, making the Danish clinical cohort an ideal testing ground for haplogroup-dependent effects.

**Predictions.** If the sink model operates in humans, haplogroup I males (more Y heterochromatin) should show weakened heterochromatin boundaries elsewhere in the genome — potentially exposing genetic variants in normally silenced regions near centromeres. Autosomal GWAS hits should differ systematically between haplogroup backgrounds, not because the autosomal variants themselves differ, but because the chromatin context determining their activity differs. This is a form of epistasis (gene-gene interaction) that standard GWAS, which ignores the Y, cannot detect.

**Preliminary evidence.** Stratifying males in the Danish clinical cohort by Y haplogroup reveals that different autosomal variants reach statistical significance depending on haplogroup background [Anonymized Reference 4]. This context-dependent association is what the sink model predicts: the same autosomal variant has different penetrance depending on whether the Y chromosome creates a strong or weak heterochromatin sink.

**An open question.** While the heterochromatin sink is established in Drosophila, whether the same mechanism operates in mammals remains debated. It is possible that mammalian Y chromosomes modulate autosomal expression through different pathways (e.g. Y-linked transcription factors). Testing this directly in clinical cohort data — rather than relying on model organisms — is a high-risk, high-reward component of this project, well suited to the exploratory spirit of the NERD programme.

### Research Plan

The project is organised around four work packages, each requiring distinct technical expertise. A key design principle is that work packages run concurrently rather than sequentially, enabling postdocs and PhD students to contribute to and benefit from results across packages in real time. Three postdocs and two PhD students will be recruited, with staggered start dates to ensure continuous supervision and cross-project synergy. Access to the large Danish ASD clinical cohort (more than 17,000 male and 5,500 female cases) is provided through an established and ongoing collaboration with a psychiatric genetics group that maintains the genotyping infrastructure and data governance. This collaboration has already produced the pilot results presented above and will continue throughout the project period, ensuring seamless data access from project start.

**Work Package 1: Genome-wide association analysis of X chromosome and Y-stratified autosomes.** This package performs the first systematic inclusion of both sex chromosomes in ASD genetic association studies. We begin with sex-stratified GWAS and X chromosome-specific association analysis (XWAS) in the Danish clinical cohort. Standard GWAS software assumes diploid autosomes; the X requires separate treatment of males (hemizygous) and females (one copy silenced by XCI). We will embed linear mixed models within population genomic models that directly estimate X-linked ancestry [Zhang et al. 2023], improving statistical power over standard approaches. Sex-stratified analyses of more than 17,000 male and 5,500 female ASD cases will quantify X-linked heritability and test whether it is concentrated in the gene sets predicted by the drive hypothesis: spermatid/neuron overlap genes, XCI escapees, and selective sweep regions.

Second, we stratify the male cohort by Y haplogroup (identified from diagnostic SNPs on the genotyping panel) and perform independent autosomal GWAS within each stratum. We will quantify missing heritability accounted for by epistatic effects of the Y chromosome using a Bayesian model that includes all haplogroups and accounts for their genealogical relationships — a hierarchical approach that borrows strength across related haplogroups rather than treating them as independent categories. We will compute the genetic correlation for ASD between haplogroup I and haplogroup R males: a low correlation would be direct evidence that the two Y backgrounds create distinct genetic architectures for the same disorder.

As results from WP3 (chromatin analysis) become available, we will associate Y-haplotype-specific GWAS hits with genes at A/B compartment borders — regions that are only accessible in either X- or Y-bearing spermatids — directly testing whether the chromatin context explains haplogroup-dependent associations. Using hitchhiking estimates from WP2, we will estimate the population-averaged polygenic risk score explained by hitchhiking-elevated frequencies of risk variants, quantifying the disease burden attributable to meiotic drive. Replication in the UK Biobank, exploiting fully sequenced X and Y chromosomes, will test generalisability. Finally, integrating results across work packages, we will develop a method to predict a male's polygenic risk score from autosomal and X-linked variants conditional on his Y chromosome repeat profile — moving from haplogroup categories to a continuous structural predictor.

**Work Package 2: Detection and characterisation of X-linked meiotic drive.** While WP1 tests statistical predictions of the drive hypothesis in clinical data, WP2 seeks direct evidence of the drive mechanism and its evolutionary history. We will run genealogy-based selection scans on the X chromosome using Relate and SINGER — tools that reconstruct ancestral recombination graphs (the genealogical history of a genomic region) from population samples and then identify loci where genealogies are shallower than expected under neutrality, indicating recent selective sweeps. Using single-cell RNA sequencing data across spermatogenesis stages, we will classify every X-linked gene by its expression trajectory: silenced at the start of meiosis and remaining off (MSCI maintained), silenced then reactivated in spermatids (MSCI escapee), or never fully silenced. The MSCI escapees — genes re-expressed after the meiotic division — are the candidate drivers. We will characterize the dual roles of each driver candidate in post-meiotic spermatogenesis and neurodevelopment using published functional data. A specific question is whether these genes encode products that remain cell-autonomous in spermatids (necessary for drive) despite the cytoplasmic bridges connecting post-meiotic cells.

The comparative dimension involves replicating selection scans in baboons using available primate genome data. Convergent positive selection on the same X-linked gene sets across primates that diverged 25–30 million years ago would be powerful evidence for intragenomic conflict rather than ecological adaptation. We will develop new genealogical statistics sensitive to transient selection episodes — the episodic flares characteristic of drive — and robust to population admixture, which confounds standard sweep detection. For each ASD-associated variant, we will estimate the change in allele frequency attributable to hitchhiking, and compute the GWAS effect size stratified by population frequency to estimate how much hitchhiking elevates the population-averaged polygenic risk score — directly quantifying the disease burden that drive imposes beyond what mutation-selection balance alone would produce.

**Work Package 3: 3D chromatin architecture in spermatogenesis.** This package provides the physical-chemistry dimension of the project, connecting the statistical genetic signals from WP1 and the evolutionary signals from WP2 to the three-dimensional organization of chromosomes in the nucleus. Hi-C (high-throughput chromosome conformation capture) measures the frequency of physical contact between every pair of genomic loci, providing a genome-wide map of chromosome folding. A fundamental feature of this organization is the division into A compartments (active, gene-rich, open chromatin) and B compartments (inactive, gene-poor, condensed). The boundaries between A and B compartments shift during spermatogenesis as chromatin is globally reorganized.

We will begin by reanalysing published single-cell Hi-C data from human spermatogenesis, mapping A/B compartment boundaries at each stage using Bayesian changepoint detection — a method that provides posterior probability distributions over boundary positions rather than sharp calls. Pilot analysis already shows that selective sweep regions on the X significantly overlap compartment boundaries in spermatogonia (p = 0.017), pachytene spermatocytes (p = 0.006), round spermatids (p = 0.005), and spermatozoa (p = 0.028) [Anonymized Reference 5]. Differentially expressed genes between X- and Y-bearing spermatids are also enriched at these boundaries (p = 0.003) [Anonymized References 5, 6].

We will identify topologically associating domains (TADs) affected at A/B borders — self-interacting chromatin neighbourhoods that constrain gene regulation — and characterize superloop structures that span compartments on the X. These large-scale chromatin loops (mediated by CTCF binding at loci such as DXZ4 and FIRRE) may enable genes in B compartments to access the transcriptional environment of A compartments, providing a structural route for MSCI escape.

Linking to the Y chromosome story: using a recently published single-cell ATAC-seq atlas of the human brain [Caglayan et al. 2024], which profiles chromatin accessibility across cell types in individuals with known Y haplogroup, we will test whether Y haplogroup influences X-linked and autosomal chromatin structure in neural tissue. If the sink model is correct, haplogroup I males should show systematically altered accessibility at euchromatin–heterochromatin boundaries. This would bridge the spermatogenesis-focused chromatin analysis with the clinical phenotype in brain.

**Work Package 4: Y chromosome repeat content and prediction.** The Y chromosome harbours large families of duplicated genes — ampliconic genes — arranged in palindromic repeat structures. These multi-copy gene families (including DAZ, CDY, RBMY, and TSPY) are essential for spermatogenesis, and their copy number varies substantially between individual men and between Y haplogroups. While WP1 uses categorical haplogroup assignments as a proxy for Y structure, this package directly characterizes the continuous variation in Y chromosome composition and develops machine learning methods to exploit it.

First, we will systematically profile ampliconic gene copy numbers and repeat content across the fully sequenced Y chromosomes now available in the UK Biobank, creating a comprehensive map of Y structural variation at population scale. This goes beyond the 43 telomere-to-telomere assemblies from Rhie et al. 2023 to thousands of individuals, enabling population-level statistical analysis.

Second, using the genealogical methods from WP2, we will characterize positive and negative selection acting on ampliconic gene counts and repeat content — testing whether certain copy-number configurations are favoured or disfavoured, and whether selection on Y repeat content connects to the meiotic drive dynamics detected on the X. If the sink model is correct, Y haplogroups with greater heterochromatin content may experience selection through their epistatic effects on autosomal gene expression — a form of indirect selection that has not previously been characterized in humans.

Third, we will develop a deep convolutional neural network (CNN) that takes a Y chromosome's repeat and ampliconic gene profile as input and predicts chromatin compartment status (A or B) at the compartment borders identified in WP3. Standard linear models cannot capture the complex, nonlinear relationship between the high-dimensional repeat landscape of the Y and chromatin state elsewhere in the genome. A CNN trained on paired Y-structure and chromatin-state data can learn these relationships end-to-end. Successful prediction would provide the strongest evidence that Y structural variation causally influences genome-wide chromatin organization, and would enable individual-level prediction of which genomic regions are epigenetically exposed — the key input to the Y-conditional risk scoring developed in WP1.

### Creativity, Ambition, and Originality

Three aspects of this project are genuinely new. First, it proposes that intragenomic conflict — a fundamental evolutionary force studied primarily in model organisms — has direct, quantifiable consequences for a specific human disorder. Meiotic drive has never been connected to human disease in this way. The prediction that selfish sex-linked genes continuously generate ASD-relevant variation represents a new explanatory framework, distinct from the standard model in which disease variants are maintained by mutation-selection balance alone.

Second, it introduces the Y chromosome as an active modulator of autosomal gene expression and ASD risk through chromatin regulation. The Y has been entirely absent from ASD genetics. If haplogroup-dependent epistasis proves real, it would provide a mechanistic basis for the male bias of ASD beyond simple X-linkage, and would explain a measurable fraction of ASD's "missing heritability" — the gap between estimated total heritability and what current GWAS variants explain.

Third, the project bridges evolutionary population genetics, comparative primate genomics, 3D chromatin biology, deep learning, and psychiatric cohort analysis — fields that rarely interact. The replication of selection scans across humans and other primates tests the generality of the conflict, while clinical cohort work tests its medical significance. The methods developed here — X-chromosome GWAS frameworks, Y-stratified heritability models, genealogical statistics for transient drive, and CNN-based chromatin prediction from Y profiles — will be applicable to any complex trait with sex bias.

### Future Applications in Life Sciences / Health Sciences

Recognizing that evolutionary conflict on the sex chromosomes shapes ASD risk could identify new risk genes and pathways invisible to current approaches. Y-haplotype-stratified polygenic risk scores may improve ASD risk prediction in males — one of the first applications of sex chromosome structural variation in clinical genomics. The CNN-based prediction framework, which maps Y chromosome structure to chromatin state, could enable stratification of ASD into subtypes with distinct aetiologies — a prerequisite for any future precision-medicine approach. The X chromosome GWAS methods developed here will be broadly useful for any GWAS — most current studies claiming "genome-wide" coverage exclude the X entirely. More broadly, the project establishes a framework for understanding how intragenomic conflict contributes to human disease, with potential relevance to other sex-biased conditions including ADHD, schizophrenia, and autoimmune disorders.

<!-- END PROJECT DESCRIPTION -->

### Abbreviations

ASD: Autism Spectrum Disorder; GWAS: Genome-Wide Association Study; XWAS: X chromosome-Wide Association Study; SNP: Single Nucleotide Polymorphism; MSCI: Meiotic Sex Chromosome Inactivation (silencing of sex chromosomes during sperm cell production); XCI: X Chromosome Inactivation (silencing of one X chromosome in female cells); PRS/PGS: Polygenic Risk Score / Polygenic Score (aggregate genetic risk from many variants); Hi-C: High-throughput Chromosome Conformation Capture (method for mapping 3D chromosome structure); ARG: Ancestral Recombination Graph (inferred genealogical history of a genomic region); HP1: Heterochromatin Protein 1 (protein that binds and maintains condensed chromatin); T2T-CHM13: Telomere-to-Telomere reference genome assembly (first complete human genome); TAD: Topologically Associating Domain (self-interacting chromatin neighbourhood); CNN: Convolutional Neural Network; LDSC-SEG: Stratified LD Score Regression (method for partitioning heritability by genomic annotation); iPSYCH: Integrative Psychiatric Research consortium (Danish clinical cohort); scATAC-seq: single-cell Assay for Transposase-Accessible Chromatin with sequencing

------------------------------------------------------------------------

## 4. LAY PROJECT DESCRIPTION (max 1,000 characters)

Autism is four times more common in boys than girls, yet no genetic explanation exists for this sex difference. We propose that evolutionary competition between the sex chromosomes contributes to autism risk. On the X chromosome, gene variants that help sperm compete spread through populations even when they slightly impair brain development — an evolutionary conflict between reproduction and neurodevelopment. On the Y chromosome, structural differences between male lineages may alter which autism risk genes are switched on. By combining evolutionary genetics, large clinical datasets with more than 22,000 autism cases, and primate genome comparisons, this project aims to uncover how sex chromosome conflict shapes the genetic basis of autism and its male bias.

------------------------------------------------------------------------

## 5. Anonymized REFERENCES (to be mapped)

- [Anonymized Reference 1] = Skov et al. 2023, "Extraordinary selection on the human X chromosome associated with archaic admixture," Cell Genomics
- [Anonymized Reference 2] = Unpublished MAGMA pilot on iPSYCH (Lundbeck-funded postdoc work)
- [Anonymized Reference 3] = Unpublished Relate analysis on 1000 Genomes (intern thesis)
- [Anonymized Reference 4] = Unpublished Y haplogroup-stratified GWAS on iPSYCH
- [Anonymized Reference 5] = Unpublished Hi-C analysis of selective sweep regions and A/B compartment borders (master thesis)
- [Anonymized Reference 6] = Unpublished single-cell RNA analysis showing SFARI ASD gene enrichment among differentially expressed genes in X- vs Y-bearing spermatids

------------------------------------------------------------------------

## NOTES FOR DISCUSSION

1. **v6 changes from v5:** Research Plan rewritten to match restructured 4-WP plan from WP_detailed_plans.md. Key changes: (a) WP1 now explicitly follows the 8-item plan with Bayesian haplogroup model, genetic correlation, and cross-WP links to A/B compartments and hitchhiking PRS; (b) WP2 adds SINGER alongside Relate, and includes explicit computation of PRS inflation from hitchhiking; (c) WP3 specifies single-cell Hi-C reanalysis and adds scATAC-seq brain atlas reference [Caglayan et al. 2024]; (d) WP4 completely rewritten — now focused on Y chromosome ampliconic gene profiling, selection on repeat content, and deep CNN for TAD chromatin prediction from Y profiles (replaces old evolutionary modelling + clinical prediction package). Formal population genetic modelling is now distributed across WP1 (Bayesian heritability, PRS methods) and WP2 (genealogical statistics, hitchhiking quantification). Added CNN, TAD, scATAC-seq, XWAS to abbreviations.

2. **Audience calibration (carried from v5):** Written for Committee for Natural and Technical Sciences — physicists, chemists, mathematicians, engineers. Every biological term defined on first use; methodology framed in terms of statistical/mathematical approaches.

3. **Anonymity risks:** The Skov et al. 2023 paper is co-authored by the applicant. The specific combination of meiotic drive + X chromosome + autism + Danish cohort could be identifying. Anonymized references must be handled carefully.

4. **Balance of risk:** WP1 and WP2 are well-supported by converging pilot data. The Y chromosome story (WP1 Y-stratification, WP3 pericentromeric boundary enrichment, WP4 CNN prediction) is higher-risk but fits the NERD "exploratory" ethos.
