# NERD 2026 Application — Working Draft v5

**Status:** Rewritten for Committee for Natural and Technical Sciences audience (physics, chemistry, mathematics, engineering)
**Deadline:** 19 February 2026, 14:00 CET

## Timeline

| Work package / Year                                      | 2027 | 2028 | 2029 | 2030 | 2031 | 2032 |
|:---------------------------------------------------------|:----:|:----:|:----:|:----:|:----:|:----:|
| WP1 GWAS of X chromosome and Y-stratified autosomes      |  X   |  X   |  X   |  X   |  -   |  -   |
| WP2 Detection and characterisation of X-linked drive      |  X   |  X   |  X   |  X   |  -   |  -   |
| WP3 3D chromatin architecture in spermatogenesis          |  -   |  X   |  X   |  X   |  -   |  -   |
| WP4 Evolutionary models and Y-based risk prediction       |  -   |  -   |  X   |  X   |  X   |  X   |
|                                                           |      |      |      |      |      |      |
| **Personnel**                                             |      |      |      |      |      |      |
| Postdoc 1 (WP1, WP2)                                     |  X   |  X   |  X   |  -   |  -   |  -   |
| Postdoc 2 (WP1, WP3)                                     |  ½   |  X   |  X   |  ½   |  -   |  -   |
| Postdoc 3 (WP2, WP3)                                     |  -   |  X   |  X   |  X   |  -   |  -   |
| PhD 1 (WP1, WP4) — 4 yr                                  |  -   |  ½   |  X   |  X   |  X   |  ½   |
| PhD 2 (WP4) — 3 yr                                       |  -   |  -   |  ½   |  X   |  X   |  ½   |

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

<!-- KM: The reader may wonder why drive has not been reported in humans and that it seems at odds the lack of evidence of transmission distortion and sex-of-offspring bias in humans: Drive is likely transient and that each driver may only produce a distortion in competitions with particular Y chromosomes and on particular autosomal backgrounds. So while each driver may strongly distort on some backgrounds its average distortion across backgrounds is often will quickly diminish or or disappear as the population allele frequencies will calibrate in favour of an average even transmission. However, before it does, the driver itself will have reached high frequency along with its hitchhiking surroundings, with the associated collateral effects. -->

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

The project is organised around four work packages, each requiring distinct technical expertise. A key design principle is that work packages run concurrently rather than sequentially, enabling postdocs and PhD students to contribute to and benefit from results across packages in real time. Three postdocs and two PhD students will be recruited, with staggered start dates to ensure continuous supervision and cross-project synergy.

**Work Package 1: Genome-wide association analysis of X chromosome and Y-stratified autosomes.** This package performs the first systematic inclusion of both sex chromosomes in ASD genetic association studies, using two complementary approaches. First, we will develop and apply X-chromosome-specific GWAS methods. Standard GWAS software assumes diploid autosomes; the X requires separate treatment of males (hemizygous) and females (one copy silenced by XCI). We will embed linear mixed models within population genomic models that directly estimate X-linked ancestry [Zhang et al. 2023], improving statistical power over standard approaches. Sex-stratified analyses of more than 17,000 male and 5,500 female ASD cases from a large Danish clinical cohort will quantify X-linked heritability and test whether it is concentrated in the gene sets predicted by the drive hypothesis: spermatid/neuron overlap genes, XCI escapees, and selective sweep regions.

Second, we will stratify the male cohort by Y haplogroup (I vs R, and subgroups) and perform independent autosomal GWAS within each stratum. Formal SNP × haplogroup interaction terms will be tested genome-wide to identify variants whose effects depend on Y background. We will compute the genetic correlation for ASD between haplogroup I and haplogroup R males — a low correlation would be direct evidence that the two Y backgrounds create distinct genetic architectures for the same disorder. ASD heritability will be estimated within each haplogroup using a Bayesian model that accounts for the genealogical relationships among haplogroups; the key test is whether heritability recovered from stratified analysis exceeds that from a pooled analysis, which would quantify how much "missing heritability" is attributable to Y-mediated epistasis.

To test the hitchhiking prediction, we will use ancestral recombination graph (ARG) methods (Relate, CLUES) to estimate allele-specific selection coefficients on the X, then evaluate the correlation between selection strength and GWAS effect size. A positive correlation — stronger selection near larger-effect ASD risk variants — would confirm that selective sweeps by driven alleles have inflated the frequencies of linked risk variants.

Replication in the UK Biobank, using fully sequenced Y chromosomes and ASD-related behavioural traits, will test generalisability. Finally, we will develop a method to predict a male's polygenic risk score from autosomal and X-linked variants conditional on his Y chromosome repeat profile — moving from haplogroup categories to a continuous structural predictor.

**Work Package 2: Detection and characterisation of X-linked meiotic drive.** While WP1 tests statistical predictions of the drive hypothesis in clinical data, WP2 seeks direct evidence of the drive mechanism and its evolutionary history. Using single-cell RNA sequencing data across spermatogenesis stages, we will classify every X-linked gene by its expression trajectory: silenced at the start of meiosis and remaining off (MSCI maintained), silenced then reactivated in spermatids (MSCI escapee), or never fully silenced. The MSCI escapees — genes re-expressed after the meiotic division — are the candidate drivers. We will test whether this set is enriched for ASD genes, for genes under positive selection, and for genes at chromatin compartment borders.

For the five triple-overlap candidate driver genes (CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1), we will characterize their dual roles in spermatogenesis and neurodevelopment using published functional data. A specific question is whether these genes encode products that remain cell-autonomous in spermatids (necessary for drive) despite the cytoplasmic bridges connecting post-meiotic cells. Proteins or mRNAs that resist sharing through bridges are the strongest driver candidates.

The evolutionary dimension of this package involves replicating selection scans across primate species. Using publicly available great ape genomes and baboon genomes accessible through consortium participation, we will run genealogy-based selection analysis (Relate/SINGER) on the X chromosome independently in each species. Convergent positive selection on the same X-linked gene sets across primates that diverged 25–30 million years ago would be powerful evidence for intragenomic conflict rather than ecological adaptation — the ecology changed, but the conflict between X and Y persists. We will develop new genealogical statistics sensitive to the transient selection episodes characteristic of drive, and robust to population admixture. For each ASD-associated variant, we will estimate how much its population frequency has been elevated by hitchhiking, and compute how much these frequency elevations inflate the population-averaged polygenic risk score — directly quantifying the disease burden attributable to drive.

**Work Package 3: 3D chromatin architecture in spermatogenesis.** This package provides the physical-chemistry dimension of the project, connecting the statistical genetic signals from WP1 and the evolutionary signals from WP2 to the three-dimensional organization of chromosomes in the nucleus. Hi-C (high-throughput chromosome conformation capture) measures the frequency of physical contact between every pair of genomic loci, providing a genome-wide map of chromosome folding. A fundamental feature of this organization is the division into A compartments (active, gene-rich, open chromatin) and B compartments (inactive, gene-poor, condensed). The boundaries between A and B compartments shift during spermatogenesis as chromatin is globally reorganized.

Using published Hi-C data from macaque spermatogenesis stages and human spermatozoa, we will map A/B compartment boundaries at each stage using Bayesian changepoint detection — a method that provides posterior probability distributions over boundary positions rather than sharp calls. Pilot analysis already shows that selective sweep regions on the X significantly overlap compartment boundaries in spermatogonia (p = 0.017), pachytene spermatocytes (p = 0.006), round spermatids (p = 0.005), and spermatozoa (p = 0.028) [Anonymized Reference 5]. Differentially expressed genes between X- and Y-bearing spermatids are also enriched at these boundaries (p = 0.003) [Anonymized References 5, 6].

We will classify boundaries as constitutive (stable across tissues) or spermatogenesis-specific (present only during meiosis), and test the prediction that meiotic drive candidates cluster specifically at dynamic boundaries — the boundaries that change during meiosis and may provide the chromatin context for MSCI escape. Extension to the Y-chromosome story links naturally: if the sink model is correct, we expect haplogroup-specific GWAS hits (from WP1) to concentrate near euchromatin–heterochromatin boundaries in brain tissue. Using the T2T-CHM13 telomere-to-telomere genome assembly — which, for the first time, resolves centromeric and pericentromeric regions completely — we will test whether haplogroup I-specific associations are enriched near pericentromeric boundaries while haplogroup R-specific ones are not.

**Work Package 4: Evolutionary modelling and Y-haplotype-based risk prediction.** This final package brings the results together through formal theory and clinical prediction. We will develop population genetic models of X-linked meiotic drive with pleiotropic effects on neurodevelopment, parameterised by empirical estimates from WP1 and WP2: the transmission distortion ratio (how strongly a driver biases segregation), the fitness cost (how much it impairs brain function), and the rate at which linked variants hitchhike. These models will predict the steady-state disease burden maintained by the evolutionary "pump" — the equilibrium frequency of ASD-relevant variation that the conflict sustains despite negative selection on the phenotype. Multi-locus extensions will model the combined effect of approximately 35 independently driven X-linked loci, predicting the expected polygenic architecture and comparing it with observed GWAS results.

The most translatable component is building a framework to predict disease variant exposure from Y chromosome structure. If the sink model holds, each Y haplotype creates a distinct chromatin landscape that determines which normally silenced risk variants become active. A male's Y repeat composition would then predict which set of risk variants is relevant to him. We will build this predictive map using iPSYCH data and validate it in the UK Biobank. The prediction model will be tested for its ability to stratify ASD cases into subtypes with distinct genetic profiles, validated against phenotypic data (age of diagnosis, co-occurring conditions, severity measures). While speculative, such stratification could eventually enable subtype-specific clinical approaches — a long-term application that motivates the fundamental work but does not constrain it.

### Creativity, Ambition, and Originality

Three aspects of this project are genuinely new. First, it proposes that intragenomic conflict — a fundamental evolutionary force studied primarily in model organisms — has direct, quantifiable consequences for a specific human disorder. Meiotic drive has never been connected to human disease in this way. The prediction that selfish sex-linked genes continuously generate ASD-relevant variation represents a new explanatory framework, distinct from the standard model in which disease variants are maintained by mutation-selection balance alone.

Second, it introduces the Y chromosome as an active modulator of autosomal gene expression and ASD risk through chromatin regulation. The Y has been entirely absent from ASD genetics. If haplogroup-dependent epistasis proves real, it would provide a mechanistic basis for the male bias of ASD beyond simple X-linkage, and would explain a measurable fraction of ASD's "missing heritability" — the gap between estimated total heritability and what current GWAS variants explain.

Third, the project bridges evolutionary population genetics, comparative primate genomics, 3D chromatin biology, and psychiatric cohort analysis — fields that rarely interact. The replication of selection scans across humans and other primates tests the generality of the conflict, while clinical cohort work tests its medical significance. The methods developed here — X-chromosome GWAS frameworks, Y-stratified heritability models, genealogical statistics for transient drive — will be applicable to any complex trait with sex bias.

### Future Applications in Life Sciences / Health Sciences

Recognizing that evolutionary conflict on the sex chromosomes shapes ASD risk could identify new risk genes and pathways invisible to current approaches. Y-haplotype-stratified polygenic risk scores may improve ASD risk prediction in males — one of the first applications of sex chromosome structural variation in clinical genomics. The variant-exposure prediction framework could enable stratification of ASD into subtypes with distinct aetiologies, a prerequisite for any future precision-medicine approach. The X chromosome GWAS methods developed here will be broadly useful for any GWAS — most current studies claiming "genome-wide" coverage exclude the X entirely. More broadly, the project establishes a framework for understanding how intragenomic conflict contributes to human disease, with potential relevance to other sex-biased conditions including ADHD, schizophrenia, and autoimmune disorders.

<!-- END PROJECT DESCRIPTION -->

### Abbreviations

ASD: Autism Spectrum Disorder; GWAS: Genome-Wide Association Study; SNP: Single Nucleotide Polymorphism; MSCI: Meiotic Sex Chromosome Inactivation (silencing of sex chromosomes during sperm cell production); XCI: X Chromosome Inactivation (silencing of one X chromosome in female cells); PRS/PGS: Polygenic Risk Score / Polygenic Score (aggregate genetic risk from many variants); Hi-C: High-throughput Chromosome Conformation Capture (method for mapping 3D chromosome structure); ARG: Ancestral Recombination Graph (inferred genealogical history of a genomic region); HP1: Heterochromatin Protein 1 (protein that binds and maintains condensed chromatin); T2T-CHM13: Telomere-to-Telomere reference genome assembly (first complete human genome); LDSC-SEG: Stratified LD Score Regression (method for partitioning heritability by genomic annotation); iPSYCH: Integrative Psychiatric Research consortium (Danish clinical cohort)

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

1. **Audience calibration (v5):** Rewritten for Committee for Natural and Technical Sciences — physicists, chemists, mathematicians, engineers. Key changes from v4: every biological term defined on first use; conceptual explanations before gene names; WP structure simplified from 9 sub-packages to 4 thematic blocks; methodology framed in terms of statistical/mathematical approaches; reduced jargon density.

2. **WP restructuring:** Following KM's outline reorganisation, the nine sub-packages have been consolidated into four thematic work packages: (1) GWAS, (2) X-linked drive detection, (3) chromatin architecture, (4) modelling and prediction. This makes the project structure more legible for non-specialists and reflects the actual expertise domains of the recruited personnel.

3. **Anonymity risks:** The Skov et al. 2023 paper is co-authored by the applicant. The specific combination of meiotic drive + X chromosome + autism + Danish cohort could be identifying. Anonymized references must be handled carefully.

4. **Balance of risk:** WP1 and WP2 are well-supported by converging pilot data. The Y chromosome story (WP1 Y-stratification, WP3 pericentromeric boundary enrichment, WP4 prediction) is higher-risk but fits the NERD "exploratory" ethos.
