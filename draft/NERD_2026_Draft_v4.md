# NERD 2026 Application — Working Draft v4

**Status:** Prose draft integrating all KM comments from v3
**Deadline:** 19 February 2026, 14:00 CET

## Timeline

| Work package / Year                                      | 2027 | 2028 | 2029 | 2030 | 2031 | 2032 |
|:---------------------------------------------------------|:----:|:----:|:----:|:----:|:----:|:----:|
| WP1.1 X chromosome GWAS                                  |  X   |  X   |  X   |  -   |  -   |  -   |
| WP1.2 MSCI escapees, drive detection, comparative        |  X   |  X   |  X   |  X   |  -   |  -   |
| WP1.3 Autosomal conflict and compensatory evolution       |  -   |  X   |  X   |  X   |  -   |  -   |
| WP1.4 3D chromatin architecture and meiotic drive         |  -   |  X   |  X   |  X   |  -   |  -   |
| WP2.1 Y haplogroup-stratified GWAS and A/B mapping        |  X   |  X   |  X   |  -   |  -   |  -   |
| WP2.2 Chromatin landscape analysis                        |  -   |  X   |  X   |  X   |  -   |  -   |
| WP2.3 Haplogroup-specific pathway and network analysis    |  -   |  -   |  X   |  X   |  X   |  -   |
| WP3.1 Population genetic models of the conflict           |  -   |  -   |  X   |  X   |  X   |  X   |
| WP3.2 Y-haplotype-based disease variant prediction        |  -   |  -   |  -   |  X   |  X   |  X   |
|                                                           |      |      |      |      |      |      |
| **Personnel**                                             |      |      |      |      |      |      |
| Postdoc 1 (WP1.1, WP1.2)                                 |  X   |  X   |  X   |  -   |  -   |  -   |
| Postdoc 2 (WP2.1, WP2.2)                                 |  ½   |  X   |  X   |  ½   |  -   |  -   |
| Postdoc 3 (WP1.3, WP1.4, WP3.1)                          |  -   |  X   |  X   |  X   |  -   |  -   |
| PhD 1 (WP2.3, WP3.1) — 4 yr                              |  -   |  ½   |  X   |  X   |  X   |  ½   |
| PhD 2 (WP3.2) — 3 yr                                     |  -   |  -   |  ½   |  X   |  X   |  ½   |

------------------------------------------------------------------------

## Potential Figures

1. Three way Venn diagram:
   - Genes subject to recent selection
   - Genes differentially expressed in X an Y spermatids
   - Genes linked to ASD
2. A/B eigenvector plot showing X chromosome A/B compartments in sperm aligned with chromosome ideogram overlaid with gene positions for genes with recent positive selection
3. Expression trajectory plots for top candidates
4. Project timeline (Gant chart)



## 1. PROJECT TITLE (max 150 characters)

**Option A (137 chars):** Selfish sex chromosomes and the genetic basis of autism: how intragenomic conflict shapes human neurodevelopment

**Option B (148 chars):** Sex chromosome conflict and autism: how meiotic drive on X and chromatin regulation by Y shape the genetic basis of neurodevelopment

**Option C (126 chars):** The evolutionary genetics of autism: selfish sex chromosomes and their consequences for human neurodevelopment

------------------------------------------------------------------------

## 2. BRIEF PROJECT DESCRIPTION (max 2,000 characters)

Autism spectrum disorder (ASD) affects nearly 3% of children in Denmark and shows a striking 4:1 male-to-female ratio that remains unexplained. This project proposes that both sex chromosomes actively contribute to ASD risk through distinct evolutionary mechanisms — the X chromosome through meiotic drive and the Y chromosome through heterochromatin-mediated gene regulation — and that understanding these mechanisms will transform our knowledge of ASD's genetic architecture.

On the X chromosome, selfish gene variants escape meiotic sex chromosome inactivation during spermatogenesis, gaining a transmission advantage even when the same variants are mildly deleterious for brain development. Because one in four neuron genes is also expressed in spermatids, this intragenomic conflict continuously pumps ASD-relevant variation into human populations. On the Y chromosome, haplogroups carrying different amounts of heterochromatin may differentially sequester chromatin-modifying proteins, creating haplogroup-dependent epistatic landscapes that modulate ASD penetrance.

The project will (1) develop population genomic methods to detect ongoing meiotic drive in humans, (2) perform the first systematic X chromosome association study of ASD using large clinical cohorts, (3) test whether Y chromosome haplogroups modulate autosomal and X-linked ASD risk, and (4) build evolutionary models integrating sex chromosome conflict with neurodevelopmental disease. The project combines methodology from evolutionary genomics with clinical cohort data, addressing questions that fall between traditional disciplines.

------------------------------------------------------------------------

## 3. PROJECT DESCRIPTION (max 30,000 characters)

<!-- BEGIN PROJECT DESCRIPTION — prose below this line counts toward the 30,000-character limit -->

### Introduction and Big Question

The 4:1 male excess in autism spectrum disorder (ASD) is one of the most robust findings in psychiatric genetics, yet its explanation remains elusive. Despite the X chromosome explaining 20% of neuro-anatomical variation relevant to ASD while constituting only 5% of the genome [Mallard et al. 2021], most genome-wide association studies have excluded the X chromosome from analysis due to its complex inheritance pattern. The Y chromosome has been ignored entirely. This project aims to fill that gap by investigating how evolutionary forces acting on both sex chromosomes shape ASD risk.

Sex chromosomes are a unique arena for evolutionary arms races because X and Y compete for transmission. An X-linked gene variant can promote its own transmission by killing Y-bearing sperm cells or by escaping meiotic sex chromosome inactivation (MSCI) — both forms of meiotic drive. Such selfish behaviour has been documented across vertebrates, insects, fungi, and plants [Zanders & Unckless 2019], and recent evidence suggests it operates in humans as well [Anonymized Reference 1]. Crucially, because one in four neuron genes is also expressed in spermatids [Matos et al. 2021], any transmission advantage conferred by spermatid expression comes at the potential cost of perturbing brain development. Evidence from a macaque half-brother study further demonstrates that paternal genotype makes an outsized contribution to heritability of a scored sub-sociality trait — a proxy for ASD that is reproducible in human children — pointing directly at spermatogenesis as the arena in which heritable neurodevelopmental risk is shaped.

On the Y chromosome, a separate mechanism operates. Haplogroups carrying different amounts of heterochromatin may differentially sequester chromatin-modifying proteins such as HP1 and SUV39H, creating haplogroup-dependent epistatic landscapes that modulate which autosomal and X-linked ASD risk variants are expressed. This "heterochromatin sink" model is well-established in Drosophila and may operate in humans, where Y chromosomal heterochromatin content varies more than tenfold between haplogroups.

This project will test the hypothesis that intragenomic conflict on the sex chromosomes — meiotic drive on the X and heterochromatin-mediated regulation by the Y — is a major unrecognized contributor to the genetic architecture of ASD and its male bias. It draws on large clinical ASD cohorts, including a Danish cohort with more than 17,000 male and 5,500 female cases plus UK Biobank data, full primate genomes that are publicly available, and baboon genomes accessible through consortium participation. The combination of evolutionary genomic methods, clinical cohort data, and comparative primate genomics positions this project to address questions that fall between traditional disciplines — exactly the kind of inter-disciplinary science the NERD programme is designed to support.

### Background — The X Chromosome and Meiotic Drive

**Recurrent selective sweeps on the human X.** The human X chromosome shows a unique pattern of recurrent bursts of strong natural selection, most recently at 45,000–55,000 years before present, where strong selection at fourteen loci dramatically affected allele frequencies for 10% of the chromosome [Anonymized Reference 1]. The scale of this phenomenon is exclusive to the X and strongly suggests selection on gene function in male spermatogenesis.

**The neuron/spermatid overlap.** A proteomic comparison shows that 80% of brain proteins are also found in testis [Matos et al. 2021]. More specifically, 24–25% of neuron genes are also expressed in male spermatids. This overlap creates the conditions for an evolutionary conflict: variants selected for their role in spermatogenesis may have collateral effects on brain development.

**The dosage-conflict hypothesis.** MSCI and somatic X chromosome inactivation (XCI) in females likely evolved from the same highly conserved mechanism of un-synapsed chromatin silencing [Hornecker et al. 2007]. Supporting this, super-loops established by CTCF binding at the XIST, FIRRE, and DXZ4 loci are crucial to proper X inactivation, and DXZ4 has been shown to function as a speciation gene in cats through its effect on MSCI — demonstrating that the same loci control both silencing mechanisms. If selfish genes escape MSCI in males, they are likely to also escape XCI in females. Because the X spends two-thirds of its time in females and one-third in males, selection will favour resolution through reduced female expression, leading to dosage deficiency in males — a direct path to male-specific neurodevelopmental risk.

**The hitchhiking hypothesis.** ASD risk variants near genes under strong meiotic drive are expected to "hitchhike" to higher population frequency than equivalent variants on autosomes. On a non-recombining chromosome, an allele linked to a driver increases in frequency even if it is mildly deleterious. Because recombination on the X is restricted to females, linked selection is particularly effective — compounded by the X chromosome's smaller effective population size (three-quarters of autosomes). The model predicts a positive correlation between selection strength and ASD risk-allele frequency, a prediction directly testable in clinical cohort data.

**Convergent evidence from comparative genomics.** The tubulin-modifying enzyme TTLL11 has been independently inactivated in both humans and gorillas but remains functional in chimpanzees. TTLL11 polyglutamylates tubulin and is required for directional sperm motility; knockout mice produce sperm that swim in circles. One hypothesis is that TTLL11 inactivation evolved as a defence against an X-linked meiotic driver that co-opted TTLL11-dependent microtubule transport across spermatid bridges to create segregation distortion. This convergent inactivation in two great ape lineages — both with low sperm competition — illustrates how meiotic drive can leave deep evolutionary signatures.

**Preliminary evidence.** Pilot analyses using MAGMA on a large clinical cohort yield results consistent with each prediction of the dosage-conflict hypothesis. In males, neuron genes co-expressed in spermatids are significantly enriched for ASD association (p = 0.030, β = 0.14), whereas neuron genes without spermatid expression show no enrichment — confirming that specifically the dual-expressed genes carry risk. ASD association in males is also elevated among genes escaping female XCI (p = 0.027, β = 0.03), consistent with MSCI-escape genes also escaping XCI. In females, ASD association is concentrated among gametologs — candidate meiotic driver genes with a Y chromosome homolog (p = 0.026, β = 0.45) [Anonymized Reference 2]. Gene Ontology analysis reveals that the neuron+spermatid overlap genes are enriched for cytosol genes (44% enrichment) and extracellular exosome genes (2-fold), while chromatin regulators are conspicuously absent — suggesting the conflict acts primarily through cytoplasmic transport and signalling rather than transcriptional regulation.

An independent analysis using the Relate genealogy method on 1000 Genomes data identified 35 positively selected genes on the X chromosome, of which 9 are ASD-associated (Fisher's exact test, p = 5.95 × 10⁻⁵). Positive selection was significantly more common in genes active in both spermatids and brain compared to brain-only genes (p = 0.002) [Anonymized Reference 3].

**Stage-resolved expression identifies MSCI escapees.** Single-cell RNA expression data through spermatogenesis allows identification of X-linked genes with post-meiotic expression — genes that must have escaped or reactivated after MSCI. These MSCI escapees are the direct candidates for meiotic drive, because post-meiotic expression is required for a gene to influence its own transmission.

**ASD genes are differentially expressed between X- and Y-bearing spermatids.** Known ASD genes from the SFARI database are significantly enriched among genes differentially expressed between X-carrying and Y-carrying spermatids (p = 0.004) [Anonymized Reference 6]. Post-meiotic spermatids share cytoplasm through intercellular bridges, meaning most transcripts are distributed equally between cell types. Genes that nonetheless show differential expression must act cell-autonomously or be regulated by the sex chromosome complement itself. The enrichment of ASD genes among these differentially expressed genes provides direct evidence that the molecular distinction between X- and Y-bearing spermatids — the distinction on which meiotic drive operates — disproportionately involves genes relevant to neurodevelopment. Furthermore, X chromosome genes differentially expressed between X- and Y-bearing spermatids are significantly enriched at A/B compartment borders (p = 0.003) [Anonymized References 5, 6], connecting the expression evidence directly to 3D chromatin architecture.

**Convergence on five candidate meiotic drivers.** Intersecting three independent gene sets — SFARI ASD genes, genes under recent positive selection (Relate), and genes differentially expressed between X- and Y-bearing spermatids — identifies five genes in the triple overlap: CDKL5, CLCN4, HUWE1, IL1RAPL1, and PTCHD1 [Anonymized References 3, 6]. Each is a well-established X-linked ASD gene. CDKL5 encodes a kinase essential for neuronal synapse formation; HUWE1 is an E3 ubiquitin ligase critical for brain development; IL1RAPL1 regulates synapse formation; PTCHD1 is among the best-replicated X-linked ASD risk genes; CLCN4 encodes a chloride channel implicated in neurodevelopmental disorders. That these five genes are simultaneously under positive selection, differentially expressed between X- and Y-bearing spermatids, and associated with ASD is precisely what the meiotic drive hypothesis predicts.

**3D chromatin architecture links selection to spermatogenesis.** Hi-C data from macaque spermatogenesis stages reveals that selective sweep regions significantly overlap A/B compartment borders — the boundaries between active and inactive chromatin domains — specifically in spermatogonia (p = 0.017), pachytene spermatocytes (p = 0.006), round spermatids (p = 0.005), and in human spermatozoa (p = 0.028) [Anonymized Reference 5]. This overlap is consistent with the hypothesis that selfish X-linked genes escape meiotic silencing by exploiting the chromatin reorganization that occurs at compartment boundaries during spermatogenesis.

### Background — The Y Chromosome and Chromatin Regulation

**The heterochromatin sink model.** The Y chromosome is the most intensively studied heterochromatin sink. In Drosophila, adding or removing Y chromosomes causes measurable genome-wide redistribution of H3K9me2/3 at pericentromeric regions [Brown et al. 2020]. The mechanism involves mass-action competition: key chromatin proteins (HP1, SUV39H) exist in limiting quantities and are titrated by heterochromatic DNA.

**Human Y haplogroup variation.** Recent assembly of 43 human Y chromosomes [Rhie et al. 2023] reveals substantial structural differences between haplogroups. Centromeric α-satellite arrays show more than 2-fold size differences (haplogroup R1b: mean 341 kbp; haplogroup I: mean 787 kbp). The DYZ1 satellite repeat varies over an order of magnitude (7–98 Mbp) and correlates with haplogroup membership.

**Chromatin regulation in ASD.** Chromatin remodelling has emerged as a central pathway in ASD genetics. High-confidence ASD risk genes include CHD8, KDM5C, SETD5, ADNP, CHD2, POGZ, and KMT5B [De Rubeis et al. 2014]. Disruption of the H3K4me3 landscape has been directly observed in autism frontal cortex samples. If Y haplogroup heterochromatin content modulates the nuclear availability of these chromatin-modifying proteins, it could create haplogroup-dependent epistatic landscapes for ASD risk.

**Predictions for haplogroup I vs R.** Haplogroup I, with its higher heterochromatin content, is predicted to cause greater sequestration of HP1, SUV39H, and H3K9 methyltransferase, leading to reduced autosomal heterochromatin formation and increased chromatin accessibility at normally silenced regions. Haplogroup R, with lower heterochromatin, should maintain stronger autosomal heterochromatin and more defined euchromatin–heterochromatin boundaries. Most Danish men carry a Y chromosome from either haplogroup I or R, making a Danish clinical cohort an ideal testing ground.

**Preliminary evidence.** Stratifying males in a large clinical ASD cohort by Y chromosome haplogroup (I vs R) reveals that different autosomal SNPs reach significance depending on haplogroup background [Anonymized Reference 4]. This suggests context-dependent epistasis mediated by Y chromosome chromatin state — precisely what the sink model predicts.

**Caveat and opportunity.** While the heterochromatin sink is well-established in Drosophila, its relevance in mammals is debated. Whether the same mechanism operates in humans, or whether mammalian Y chromosomes modulate autosomal expression through different pathways (e.g. ZFX/ZFY transcription factors), remains an open question. Testing this directly in human ASD cohorts represents a high-risk, high-reward component — exactly the kind of exploratory science the NERD programme is designed to support.

### Research Plan

#### Work Package 1: The X Chromosome — Detecting Meiotic Drive and Its Consequences (Years 1–4)

**WP1.1: X chromosome GWAS in clinical cohorts.** The first aim is to perform the first large-scale, methodologically rigorous X chromosome association study for ASD. We will develop population genomic methods sensitive to rare X-linked variants by embedding standard linear mixed models for GWAS in population genomic models of cohort ancestry [Zhang et al. 2023], directly modelling ancestry along the X chromosome to increase power over standard imputation. We will run sex-stratified GWAS on a large clinical cohort with more than 17,000 male cases and more than 5,500 female cases, using both dosage-based and hard-call genotypes to evaluate the impact of dosage compensation modelling. ASD heritability on the X will be estimated using stratified LD score regression (LDSC-SEG), partitioning heritability across nested gene sets: spermatid/neuron overlap, XCI escapees, gametologs, and genes in selective sweep regions. X-linked polygenic risk scores (PRS) will be built and tested for predictive power beyond autosomal PRS — quantifying what the field misses by excluding the X. The hitchhiking prediction will be tested directly: we will use RELATE and CLUES to compute the likelihood ratio of positive selection at each X-linked SNP, then evaluate the expected positive correlation between selection strength and risk-allele frequency. Data from the PGC and SPARK cohorts will be pooled to increase sample size and population diversity.

**WP1.2: Identifying MSCI escapees, detecting meiotic drive, and comparative genomics.** Using single-cell RNA expression data across spermatogenesis, we will identify X-linked genes with post-meiotic expression (MSCI escapees) and test whether they are enriched for ASD-associated genes. Expression trajectories of candidate genes through spermatogenesis stages will characterize the silencing-reactivation pattern that constitutes direct evidence of MSCI escape. The five triple-overlap genes (CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1) will be the starting point for functional characterization, investigating their involvement in selective or directional transport across cytoplasmic bridges between X- and Y-bearing spermatids, uneven distribution of mitochondria to X- and Y-carrying spermatids, and differential storage of mRNA in chromatoid bodies — the cell-biological mechanisms through which drive could operate at the molecular level. Genealogy-based methods (Relate, ARG inference) will be applied to X chromosome data from multiple population datasets to detect signatures of ongoing drive. Crucially, we will replicate selection scans in macaque and baboon genomes to test whether the same gene sets are repeatedly targeted by strong positive selection across primates. Convergent selection on the same X-linked genes across independently evolving primate lineages would provide powerful evidence that intragenomic conflict, rather than ecological adaptation, is the driving force. Full primate genomes are publicly available, and baboon genomes are accessible through consortium participation.

**WP1.3: The neuron/spermatid conflict — autosomal impact and compensatory evolution.** This work package extends the conflict hypothesis to the autosomes, testing whether autosomal genes functionally linked to X- or Y-linked genes are enriched for signatures of selection — the expected signature of continuous compensatory evolution against transmission distortion by drivers. Scans for selection on autosomes will use publicly available 1000 Genomes data, asking whether compensatory selection is detectable as enrichment of positive selection among autosomal interactors of X-linked driver candidates. We will estimate pathogenicity of segregating variants to test whether genes involved in spermatogenesis accumulate more pathogenic variants than expected, and investigate whether autosomal genes co-opted by candidate drivers are more often disrupted by in-frame stop codons. The shared gene functions and pathways between spermatids and neurons that are under selection will be characterized; preliminary GO analysis shows enrichment for cytosol and extracellular exosome genes, suggesting the conflict operates through cytoplasmic transport and signalling. The Doublecortin (DCX) family will serve as a prime case study: DCX is the master X-linked regulator of neuronal migration through microtubule stabilization and is also required for sperm motility, shows strong recent selection, and is associated with ASD. TTLL11 inactivation will be investigated as a model for organismal defence against meiotic drive. Multivariate PGS regression will test for differential genetic load across ASD subtypes, evaluating whether the meiotic drive gene set explains phenotypic variation within the autism spectrum.

**WP1.4: 3D chromatin architecture and meiotic drive.** We will analyze Hi-C data across spermatogenesis stages to characterize how A/B compartment borders relate to selfish gene regions, applying Bayesian compartment border detection to precisely map boundaries in spermatogonia, spermatocytes, spermatids, and spermatozoa. Properties of candidate genes located at fluid A/B compartment borders will be characterized, with particular attention to understanding which borders are shared across tissues and which change specifically during spermatogenesis — a distinction critical for understanding whether drive exploits constitutive or dynamic chromatin features. We will test whether genes escaping MSCI cluster at compartment boundaries, providing a mechanistic link between chromatin reorganization and meiotic drive, and extend the analysis to human spermatogenesis Hi-C data as datasets become available.

#### Work Package 2: The Y Chromosome — Haplogroup-Dependent Epistasis (Years 1–5)

**WP2.1: Y haplogroup-stratified GWAS and chromatin-state mapping.** We will assign Y haplogroup to each male in the clinical cohort using Y-chromosomal SNPs, validate assignments against the T2T-CHM13 Y chromosome assemblies, partition the cohort by haplogroup (I vs R and subgroups), and run independent autosomal association tests. Formal SNP × haplogroup interaction terms will be tested genome-wide using GEM (Gene-Environment interaction in Millions), treating haplogroup as the "environment." Beyond standard GWAS, we will identify autosomal genes whose A/B compartment status in males is determined by Y chromosome haplotype — the direct mechanistic test of the sink model. Heritability in males will be estimated using a Bayesian method that stratifies by Y haplotype, quantifying how much of the heritable component of ASD is haplogroup-dependent. The analysis will be extended to cohorts where deeply diverged Y haplogroups provide a more extreme test of the sink hypothesis.

**WP2.2: Chromatin landscape analysis.** Haplogroup-specific significant SNPs will be mapped to Roadmap Epigenomics chromatin states for brain tissues (prefrontal cortex, hippocampus, cerebellum). We will test enrichment in H3K9me3-marked heterochromatin, H3K4me3-marked promoters, and enhancers using LDSC-SEG and GoShifter. The key prediction — that SNPs near euchromatin–heterochromatin boundaries should show stronger effects in haplogroup I background — will be tested using the T2T-CHM13 assembly with its complete pericentromeric annotations. If methylation data with Y haplogroup information is available, haplogroup-associated differentially methylated regions at pericentromeric boundaries and ASD gene promoters will be tested using minfi and DMRcate.

**WP2.3: Haplogroup-specific pathway and network analysis.** Haplogroup-specific gene lists will be mapped to GO terms, KEGG pathways, and Reactome, focusing on chromatin modification, synaptic function, and neurodevelopmental pathways. We will test whether genes significant in different haplogroup backgrounds converge on the same pathways despite being different genes — evidence for the same biology being affected through different chromatin routes. Haplogroup-specific PRS will be built and tested for improved ASD risk prediction against unstratified models, directly quantifying what fraction of heritability that appears missing in a combined analysis is recovered when stratifying by Y haplogroup.

#### Work Package 3: Evolutionary Modelling and Clinical Prediction (Years 3–6)

**WP3.1: Population genetic models of the neuron/spermatid conflict.** We will develop formal population genetic models of X-linked meiotic drive with pleiotropic effects on neurodevelopment, parameterising the transmission ratio distortion, the fitness cost in brain function, and the frequency-dependent dynamics. The models will predict the steady-state ASD risk generated by the evolutionary pump: given realistic parameters for drive strength and fitness costs, what equilibrium frequency of ASD-relevant variation does the conflict maintain? The expected polygenic architecture under the conflict model will be compared with neutral expectations, predicting excess common variation at dual-expressed genes, higher risk-allele frequencies near selected loci, and enrichment of ASD association on the X. We will quantify PGS contributed by autosomal genes whose A/B compartment status in males is determined by Y haplotype, testing whether these account for the higher ASD prevalence in haplogroup I versus R males. A key empirical analysis will estimate how many causal variants on the X chromosome have hitchhiked to frequencies much higher than their detrimental effect would predict — the expected signature when meiotic drive trumps the detrimental neurodevelopmental effect of linked variants.

**WP3.2: Y-haplotype-based disease variant prediction.** The most clinically translatable aim of the project is to predict the subset of disease variants exposed in a male based on his Y chromosome haplotype and repeat composition. If the sink model is correct, a male's Y haplotype determines which chromatin boundaries are relaxed, which in turn determines which risk variants in normally silenced regions become active. This creates the possibility of stratifying ASD cases by Y-haplotype-predicted variant exposure, potentially isolating subtypes of ASD with separate aetiologies. Such stratification could enable better and earlier support for ASD children and inform development of treatments targeted to subtypes — for example, subtype-specific responses to Lamotrigine (which modulates glutamate signalling) or Oxytocin (which affects social cognition pathways). While these clinical applications remain speculative at this stage, establishing the Y-haplotype-to-variant-exposure mapping is a concrete deliverable of this project and a prerequisite for any future precision-medicine approach to sex-chromosome-mediated ASD risk.

<!-- END PROJECT DESCRIPTION -->

### Creativity, Ambition, and Originality

This project is original in three ways. First, it proposes that intragenomic conflict — a fundamental evolutionary force — has direct consequences for human disease. While meiotic drive has been studied extensively in model organisms, no one has connected it to a specific human disorder. The idea that selfish X-linked genes continuously pump ASD-relevant variation into human populations represents a paradigm shift in understanding the genetic architecture of neurodevelopmental disorders.

Second, it introduces the Y chromosome as an active modulator of ASD risk through chromatin regulation. The Y has been entirely ignored in ASD genetics. If haplogroup-dependent epistasis proves real, it would explain a substantial portion of ASD's "missing heritability" and provide a mechanistic basis for its male bias beyond simple X-linkage. The new aim of predicting exposed disease variants from Y haplotype extends this idea toward clinical utility.

Third, it bridges evolutionary genomics, comparative primate genomics, and psychiatric genetics — fields that rarely interact. The replication of selection scans across humans, macaques, and baboons tests the generality of the conflict, while the clinical cohort work tests its medical significance. The methods and conceptual framework developed here will be applicable to any complex trait with sex bias.

### Future Applications in Life Sciences / Health Sciences

The proposed research has potential future applications across several domains. By revealing how evolutionary forces on sex chromosomes shape ASD risk, this project could identify entirely new classes of risk genes and pathways that current approaches miss. Haplogroup-stratified polygenic risk scores may substantially improve ASD risk prediction in males — one of the first examples of incorporating sex chromosome variation into clinical genomics. The Y-haplotype-based variant prediction framework could enable stratification of ASD into subtypes with distinct aetiologies, opening the door to earlier diagnosis and subtype-targeted interventions. The X chromosome association methods developed here will be broadly useful, enabling systematic inclusion of X-linked markers in GWAS of any trait — addressing a major gap in a field that routinely excludes the X from studies claiming to be "genome-wide." More broadly, the project establishes a framework for understanding how intragenomic conflict contributes to human disease, applicable to other sex-biased disorders including ADHD, schizophrenia, and autoimmune conditions. Understanding the X chromosome's disproportionate contribution to ASD is crucial for improved diagnosis and support of neurodivergent individuals.

### Abbreviations

ASD: Autism Spectrum Disorder; GWAS: Genome-Wide Association Study; MSCI: Meiotic Sex Chromosome Inactivation; XCI: X Chromosome Inactivation; SNP: Single Nucleotide Polymorphism; PGS: Polygenic Score; PRS: Polygenic Risk Score; GO: Gene Ontology; HP1: Heterochromatin Protein 1; ARG: Ancestral Recombination Graph

------------------------------------------------------------------------

## 4. LAY PROJECT DESCRIPTION (max 1,000 characters)

Autism is increasingly common, affecting nearly 3% of Danish children, and is four times more prevalent in boys than girls. This sex difference remains unexplained. This project investigates whether "selfish" genes on the sex chromosomes contribute to autism risk. On the X chromosome, gene variants that enhance sperm production spread through populations even when they mildly impair brain development. On the Y chromosome, structural differences between lineages may alter which autism risk genes become active. By combining evolutionary genetics with analysis of large clinical datasets, this project aims to uncover new genetic mechanisms underlying autism and its male bias.

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

1. **Framing priority:** The project is framed as fundamental evolutionary genomics (natural science) with applications in health — not as biomedicine.

2. **Anonymity risks:** The Skov et al. 2023 paper is co-authored by the applicant. The specific combination of meiotic drive + X chromosome + autism + Danish cohort could be identifying. Anonymized references must be handled carefully.

3. **Balance of risk:** WP1 is now very well-supported by five converging lines of pilot data. WP2 is higher-risk but fits the NERD "exploratory" ethos. WP3 bridges modelling and clinical translation.

4. **Changes from v3:**
   - WP1.2 expanded with comparative primate genomics (macaque/baboon selection scans) and functional characterization of driver mechanisms (bridge transport, mitochondrial asymmetry, chromatoid bodies)
   - WP1.3 expanded with autosomal compensatory evolution tests, pathogenicity estimation, and co-option analysis
   - WP1.4 expanded with tissue-shared vs spermatogenesis-specific border characterization
   - WP2 moved up in the narrative to directly follow WP1
   - WP2.1 expanded with A/B compartment mapping by Y haplotype and Bayesian stratified heritability
   - Old WP3.2 (chromatin simulation) and WP3.3 (X-Y integration) removed as infeasible
   - New WP3.2 added: Y-haplotype-based disease variant prediction and ASD subtype isolation
   - WP3.1 expanded with PGS quantification from Y-dependent A/B genes and hitchhiking frequency analysis
   - Introduction now references macaque half-brother study, iPSYCH/UK Biobank access, and primate genome availability
   - Timeline table filled in
