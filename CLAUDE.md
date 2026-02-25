# CLAUDE.md — NERD 2026 Grant Application Project

## Project overview

Grant application by Kasper Munch (Associate Professor, Aarhus University BiRC) investigating how evolutionary forces on sex chromosomes — meiotic drive on X, heterochromatin sink on Y — shape the genetic architecture of autism spectrum disorder (ASD).

The NERD 2026 application (Novo Nordisk Foundation, up to DKK 16M over 6 years) was submitted 19 Feb 2026. The project is now being expanded for submission to other funding bodies.

## Key files

- `NERD_2026_Draft.md` — Submitted NERD application (Title, Brief Description, Project Description, Lay Description)
- `references.bib` — BibTeX file for all cited works
- `meiotic_drive_section.md` — Standalone expanded description of meiotic drive theory (for reuse in other applications)
- `research.md` — PI's research summary
- `WP_detailed_plans.md` — Detailed work package plans
- `wishlist.md` — Pilot result wishlist with status
- `TODO.md` — Task tracking (partially outdated after submission)
- `gantt.yaml` — Authoritative project timeline (WP1 2027–2030, WP2 2027.5–2029.5, WP3 2028.5–2031.5, WP4 2029.5–2033)
- `budget/NNF_NERD_budget_template_NAT_08012026_2025_project_supplement_MBG.xlsx` — Budget (DKK 15,407,395 total)
- `source-docs/` — Reference PDFs and background research
- `figures/` — Application figures (4 figures: drive model, Venn+expression, A/B compartments with RBMX2, Gantt chart)

## Core conceptual points (corrections from Kasper)

These are critical framings that must be respected in all writing:

1. **Y heterochromatin is NOT a parallel study** — it is the MEANS by which the Y exerts control over transmission. The Y's bulk of repetitive DNA is its weapon in the meiotic drive arms race. Never frame Hypothesis 2 as "in parallel" to Hypothesis 1.

2. **The project is NOT about a technical GWAS fix** — it proposes an entirely new evolutionary framework (intragenomic conflict + non-adaptive evolution) for understanding ASD. Do not frame the innovation as "including X chromosomes in GWAS."

3. **Hypothesis 2 presumes repeat content evolves as a means to control transmission** — Y chromosome repeat content is not just structural variation; it evolves under selection as the Y's weapon in the arms race.

4. **Terminology**: Use "hitchhiking genomic neighborhood" as the single lay term replacing both "selective sweep" and "hitchhiking" when writing for non-specialists. Do not use two separate popgen terms.

5. **The @lemos2008 citation** (Lemos, Araripe & Hartl 2008, Science) is cited for Y chromosome kmer profiling predicting chromatin state — but Kasper was looking for a specific Andy Clark paper on this topic. The citation may need replacing.

## Application structure (NERD 2026)

- **Title**: max 150 chars
- **Brief Description**: max 2,000 chars
- **Project Description**: max 30,000 chars (between BEGIN/END comment markers in the .md)
- **Lay Description**: max 1,000 chars
- **Literature References**: max 8,000 chars (separate field)
- **Illustrations**: max 4 (PNG/JPG/BMP, ≤1050×1650px)

Assessment criteria: scientific quality, creativity, novelty, ambition, feasibility, dedication (time).
Audience: Committee for Natural and Technical Sciences (physicists, chemists, mathematicians, engineers).
Stage 1 is anonymous.

## Working with Kasper — editing workflow

- **Go paragraph by paragraph or section by section**, proposing changes for acceptance/rejection.
- **Write proposed text to a temporary named markdown file** that Kasper can review.
- **Only move content into the main document when Kasper explicitly accepts it.**
- This applies to both modifying existing texts and creating new documents from outlines.
- Never rewrite entire documents at once.
- Kasper is selective — he may accept only some of multiple proposed changes. Always ask or present options individually when there are independent changes.

## Abbreviations used in the project

ASD: Autism Spectrum Disorder; GWAS: Genome-Wide Association Study; MSCI: Meiotic Sex Chromosome Inactivation; XCI: X Chromosome Inactivation; Hi-C: High-throughput Chromosome Conformation Capture; CNN: Convolutional Neural Network; SFARI: Simons Foundation Autism Research Initiative (gene database); PRS: Polygenic Risk Score; WP: Work Package

## Work packages (NERD version)

- **WP1**: Y-stratified GWAS including X chromosome (Danish cohort 17K male + 5.5K female cases; UK Biobank replication)
- **WP2**: Detection and characterization of meiotic drive (genealogy-based selection scans, hitchhiking quantification, baboon replication)
- **WP3**: 3D chromatin architecture in spermatogenesis (Hi-C reanalysis, Bayesian changepoint detection, scATAC-seq brain atlas)
- **WP4**: Y chromosome repeat content and prediction (k-mer profiling, Brownian motion model, deep CNN for chromatin prediction)

## Personnel

PhD1, PD1, PhD2, PD2, PD3 — staggered start dates across WPs.
