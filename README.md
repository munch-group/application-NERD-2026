# NERD 2026 Application — Selfish Sex Chromosomes and Autism

**Deadline:** 19 February 2026, 14:00 CET
**Applicant:** Kasper Munch, Lektor, BiRC, Aarhus University
**Grant:** Novo Nordisk Foundation NERD Programme — up to DKK 16M over 6 years

## For Claude: Session Briefing

This repository contains all materials for the NERD 2026 grant application. To resume work, read this README first, then `draft/NERD_2026_Draft_v1.md`, then `wishlist.md`.

## Repository Structure

```
├── README.md                  ← This file (session briefing + all context)
├── draft/
│   └── NERD_2026_Draft_v1.md  ← Working draft of the full application
├── wishlist.md                ← Prioritised pilot result wishlist
├── source-docs/
│   ├── Guidelines-NERD-Programme-2026.pdf    (NERD 2026 official guidelines)
│   ├── R436-2023-1175.pdf                    (Lundbeck application — successful)
│   ├── THESIS_ARIADNA_SAEZ.pdf               (Intern thesis — Relate analysis)
│   ├── Y_Chromosome_Haplogroup_Autism_Research_Report.md
│   └── heterochromatin_sink_review.md
├── figures/
│   └── (Kasper to add Venn diagram image here)
└── gene-lists/
    └── venn_overlaps.md       ← Complete gene lists from the Venn diagram
```

## The Grant

Novo Nordisk Foundation NERD (New Exploratory Research and Discovery) Programme 2026. Fundamental research in natural/technical sciences with applications in life/health/sustainability. **"Bioscience and basic biomedicine" is explicitly OUT of scope** — must frame as evolutionary genomics.

**Anonymous first round:** Only the Proposal tab is seen (title, brief description, project description, lay description, illustrations, literature references). No identifying information.

**Second round (if invited):** CV, publication list, summary of own research, budget.

### Character Limits

| Section | Limit |
|---------|-------|
| Project Title | 150 chars |
| Brief Project Description | 2,000 chars |
| Project Description | 30,000 chars |
| Literature References | 8,000 chars |
| Lay Project Description | 1,000 chars |
| CV | 4,000 chars |
| Publication list | 5,000 chars |
| Summary of own research | 2,000 chars |
| Illustrations | max 4, JPG/PNG/BMP, max 1050×1650 px |

### Assessment Criteria

Scientific quality, creativity, novelty, ambition, feasibility, dedication (time commitment)

---

## The Project

### Core Hypothesis
Both sex chromosomes contribute to ASD risk: the X through meiotic drive, the Y through heterochromatin-mediated chromatin regulation.

### X Chromosome — Meiotic Drive (WP1, Years 1–4)
Selfish X-linked gene variants escape meiotic sex chromosome inactivation (MSCI) during spermatogenesis, gaining a transmission advantage. Because 24-25% of neuron genes are also expressed in spermatids, variants selected for spermatogenesis can harm brain development. Three sub-hypotheses:

1. **Dosage-conflict:** MSCI and XCI share mechanisms. Genes escaping MSCI also escape XCI in females → reduced female expression → male dosage deficiency → male-specific neurodevelopmental risk
2. **Hitchhiking:** ASD risk variants near selfish genes rise to higher frequency due to linked selection on the non-recombining (in males) X chromosome
3. **Compartment border escape:** Selfish genes exploit chromatin reorganisation at A/B compartment borders during spermatogenesis to escape silencing

### Y Chromosome — Heterochromatin Sink (WP2, Years 2–5)
Y haplogroups carry different amounts of heterochromatin (haplogroup I > haplogroup R). Heterochromatin sequesters chromatin-modifying proteins (HP1, SUV39H) in limiting quantities. Different haplogroups create different epistatic landscapes for ASD risk. Well-established in Drosophila; weaker evidence in mammals (high-risk, high-reward).

### Evolutionary Modelling (WP3, Years 3–6)
Formal population genetic models of meiotic drive with pleiotropic neurodevelopmental effects. Chromatin factor titration simulations. X–Y interaction models.

---

## Preliminary Results (all confirmed)

### MAGMA pilot on iPSYCH (Anonymised Reference 2)
- Neuron+spermatid genes enriched for ASD in males: p=0.030, β=0.14
- XCI escape genes enriched for ASD in males: p=0.027, β=0.03
- Gametologs enriched for ASD in females: p=0.026, β=0.45
- GO: cytosol genes 44% enriched, exosome genes 2×, chromatin genes absent in neuron+spermatid set
- Key dual-function XCI-escaper genes: DDX3X, HUWE1, UBA1, SYAP1, FMR1

### Relate selection scan on 1000 Genomes (Anonymised Reference 3)
- 35 positively selected genes on X, 9 are ASD-associated (p=5.95×10⁻⁵)
- Selection more common in spermatid+brain vs brain-only genes (p=0.002)
- Key ASD genes under selection: CASK, HUWE1, FRMPD4, IL1RAPL1, CLCN4, CDKL5, DMD, PTCHD1

### Single-cell RNA spermatogenesis (Anonymised Reference 6)
- SFARI ASD genes enriched among DEGs in X- vs Y-bearing spermatids: **p=0.004**
- Those DEGs enriched at A/B compartment borders: **p=0.003**
- **Triple-overlap Venn diagram:** 5 genes in SFARI × selection × X/Y spermatid DEGs: **CDKL5, CLCN4, HUWE1, IL1RAPL1, PTCHD1**
- All pairwise overlaps significant (see `gene-lists/venn_overlaps.md`)

### Hi-C compartment borders (Anonymised Reference 5)
- ECH90 selective sweep regions overlap A/B borders in macaque spermatogonia (p=0.017), pachytene spermatocytes (p=0.006), round spermatids (p=0.005), human spermatozoa (p=0.028)
- Bayesian method developed for precise border detection
- Hi-C thesis: https://munch-group.github.io/hic-spermatogenesis/thesis
- Bayesian method: https://munch-group.github.io/hic-spermatogenesis

### Y haplogroup-stratified GWAS (Anonymised Reference 4)
- Stratifying iPSYCH males by Y haplogroup (I vs R) changes which autosomal SNPs reach significance
- Qualitative result confirmed; specific statistics not yet shared

---

## Anonymised Reference Mapping

| Reference | Source |
|-----------|--------|
| Anonymised Reference 1 | Skov et al. 2023, Cell Genomics — "Extraordinary selection on the human X chromosome" (co-authored by Kasper) |
| Anonymised Reference 2 | Unpublished MAGMA pilot on iPSYCH (Lundbeck-funded postdoc work by Shannon D'Urso) |
| Anonymised Reference 3 | Unpublished Relate analysis on 1000 Genomes (intern thesis by Ariadna Saez Gómez) |
| Anonymised Reference 4 | Unpublished Y haplogroup-stratified GWAS on iPSYCH |
| Anonymised Reference 5 | Unpublished Hi-C analysis of selective sweep regions and A/B compartment borders (master thesis by Søren Jørgensen) |
| Anonymised Reference 6 | Unpublished scRNA-seq analysis: SFARI enrichment among DEGs in X- vs Y-bearing spermatids |

## Key Published References

- Skov et al. 2023, Cell Genomics — "Extraordinary selection on the human X chromosome associated with archaic admixture"
- Mallard et al. 2021 — X explains 20% of neuro-anatomical variation relevant to ASD
- Matos et al. 2021 — 80% of brain proteins also found in testis
- Hornecker et al. 2007 — MSCI/XCI shared mechanism
- Zanders & Unckless 2019 — meiotic drive across taxa
- Zhang et al. 2023 — population genomic GWAS ancestry models
- Brown et al. 2020 MBE — Y heterochromatin sink in Drosophila
- Rhie et al. 2023 — 43 human Y chromosome assemblies, haplogroup structural variation
- De Rubeis et al. 2014 — chromatin remodelling genes in ASD (CHD8, KDM5C, SETD5, ADNP, CHD2, POGZ, KMT5B)

---

## Group Members (for budget planning)

**Current:**
- Kasper Munch (PI, Lektor — permanent, salary NOT coverable by NERD)
- Shannon D'Urso (postdoc, Lundbeck-funded through ~2025, works on X chromosome GWAS and meiotic drive/ASD)
- Johan Christensen Ulstrup (master student, baboon X chromosome)

**Former (relevant):**
- Erik Fogh Sørensen (PhD, X chromosome evolution and hybrid incompatibilities)
- Søren Jørgensen (master, Hi-C spermatogenesis — Anonymised Ref 5)
- Ariadna Saez Gómez (intern, Relate selection analysis — Anonymised Ref 3)
- Davide Capozzi (intern, meiotic drive dynamics)
- Tobias Røikjer (master, phase-type distributions)

**NERD budget should cover:** 1-2 PhD students, 1-2 postdocs, computing, travel. No PI salary. No overhead.

## iPSYCH Cohort Details

Shannon is working with iPSYCH2015 release:
- Males: ~17,300 cases / ~22,500 controls
- Females: ~5,600 cases / ~22,200 controls
- GRCh37 coordinates
- X chromosome imputed separately (no PARs in current imputation)
- Y haplogroups assignable from Y-chromosomal SNPs

---

## Decisions Made

- Topic: meiotic drive and autism (NOT phase-type distributions)
- Phase-type distributions: NOT included in WP3
- Framing: fundamental evolutionary genomics with health applications (NOT biomedicine)
- Three work packages: WP1 (X/meiotic drive, years 1-4), WP2 (Y/haplogroup epistasis, years 2-5), WP3 (modelling, years 3-6)

## Remaining Tasks

1. **Convert draft from outline to prose** — Project description is currently outline/bullet format. Convert to flowing prose. Currently ~24,900 chars; target close to 30,000 chars.
2. **Design 4 figures/illustrations** — Proposed:
   - Figure 1: Conceptual overview (two mechanisms converging on ASD)
   - Figure 2: Venn diagram (produced by Kasper, add to `figures/`) + expression trajectories if available
   - Figure 3: Hi-C compartment border results + Y haplogroup data
   - Figure 4: Timeline / work package Gantt chart
3. **Prepare budget breakdown** — DKK 16M over 6 years
4. **Draft CV section** (4,000 chars)
5. **Draft publication list** (5,000 chars)
6. **Draft summary of own research** (2,000 chars)
7. **Literature references** (8,000 chars) — with anonymised self-citations
8. **Final anonymity check** — no identifying information in Proposal tab

## Relevant Notion Pages

- **Research Tree:** https://www.notion.so/208fd1e7c2e180ee9aacc44071c02889
- **People:** https://www.notion.so/267fd1e7c2e1801ebf85cad3a6a26c7c
- **X Chromosome GWAS (Shannon's work log):** https://www.notion.so/2edfd1e7c2e1810996ecc4b263450ada
- **GWAS with chrY haplotypes:** https://www.notion.so/244fd1e7c2e18032819ed1eb7db4ff57
- **Research strategy:** https://www.notion.so/208fd1e7c2e180beb13ee46807c546e9
- **Pilot wishlist (Notion copy):** https://www.notion.so/306fd1e7c2e181f18ff0e7aadaf58248
