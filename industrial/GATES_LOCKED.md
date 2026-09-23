# GATES LOCKED - industrial/ slice (builder 21 replacement), project enzyme-mining-pollution
Locked 2026-09-23 18:55 IST, BEFORE any mining/results. Committed standalone per fleet rule.

Scope: industrial-pollutant-degrading enzymes mined from public metagenomes (MGnify
soil + marine-sediment catalogues) + UniProt background corpora + decoy proteomes.
Amended bar (shared with B19/B20): mining platform + validated specificity gate,
NOT another curation DB. Two mechanistically distinct classes:
- HLD: hydrolytic dehalogenases (haloalkane dehalogenase EC 3.8.1.5 + haloacid
  dehalogenase EC 3.8.1.2) - chlorinated industrial solvents/intermediates.
- RHD: Rieske non-heme aromatic ring-hydroxylating dioxygenase alpha subunits
  (EC 1.14.12.x: benzene/toluene 1.14.12.3/.11, naphthalene 1.14.12.12, biphenyl
  1.14.12.18 and related) - petroleum/PAH/BTX industrial pollutants.

## Prior-art verdict: CROWDED (specific contribution named, not a duplicate)
Closest works:
1. EnzymeMiner 2.0 (NAR 2025, https://academic.oup.com/nar/article/54/W1/W257/8675561) -
   automated sequence mining platform, generic, no frozen-threshold specificity
   gating against documented false-positive modes.
2. Metagenome-derived haloalkane dehalogenases with novel catalytic properties
   (PubMed 28674849, 2017) - functional HLD mining, small scale, no discriminator.
3. Characterization of novel PAH dioxygenases from contaminated-soil metagenomic DNA
   (AEM 2014, doi:10.1128/aem.01883-14) and similar PAH-RHD mining works.
Contribution: extending the fleet's frozen-threshold, negative-mode-gated mining
platform (validated in pesticides/) to two chemistries that are NOT hydrolase
esterases - SN2 hydrolytic dehalogenation and Rieske monooxygenation - plus the
cross-class unified vetted-candidate database (B19+B20+B21) with per-candidate
evidence rows. Not done before as a gated, byte-locked pipeline.

## Reference set (assembled + byte-locked 2026-09-23; data/reference_set.csv)
- HLD: 42 VERIFIED (33 haloalkane + 9 haloacid) + 12 THIN = 54
- RHD: 51 VERIFIED + 12 THIN = 63. Total 117 (93 VERIFIED / 24 THIN).
- Tiers per B20 convention: VERIFIED = Swiss-Prot reviewed w/ target-activity
  evidence; THIN = TrEMBL annotation-only propagation, tracked separately, never
  silently promoted. Gates measured on VERIFIED held-out members only.
- Split BEFORE mining: seeded (seed=20260923) 70/30 stratified by (class, family),
  frozen in data/split.json (82 train / 35 heldout; 31 of heldout VERIFIED).
  Held-out sequences never used for profiles.

## GATE 1 - Reference recall >=90% (parent-set)
Pipeline (family profile HMMs + motif checks) must recover >=90% of HELD-OUT
VERIFIED enzymes hidden in the mining corpus (MGnify MAGs + bg corpora + decoy
proteomes). Positive control first: 100% recovery of TRAINING VERIFIED enzymes.
Operating thresholds fixed by training-set cross-validation ONLY, before the
held-out test. No post-hoc tuning; a failed gate is reported as a failure.

## GATE 2 - Specificity discriminator vs documented false-positive modes
Labeled negatives byte-locked in data/negative_set.csv (70 total; 34 tune / 36 eval):
  (a) 20 HAD-superfamily phosphatases/nucleotidases - the classic haloacid
      dehalogenase-fold false positive (HAD-fold != halidohydrolase);
  (b) 15 bacterial epoxide hydrolases - alpha/beta-hydrolase-fold neighbors of
      haloalkane dehalogenases with no dehalogenation activity;
  (c) 20 Rieske [2Fe-2S] ferredoxins - electron-transfer components of the same
      systems, no oxygenase domain;
  (d) 15 intra/extradiol ring-CLEAVAGE dioxygenases (EC 1.13.11.x) - aromatic
      dioxygenases of different chemistry (cleavage, not ring hydroxylation).
Discriminator emits CONFIDENT / PROMISCUOUS-TIER / REJECT. At frozen thresholds:
precision >=90% on held-out labeled VERIFIED-positives-vs-negatives AND AUROC >=0.9.

## GATE 3 - Novelty
Claimed novel candidates: <80% pairwise identity to any VERIFIED reference enzyme
(3-mer Jaccard proxy, flagged) AND absent from BRENDA/literature curation for the
claimed activity (checked, not assumed).

## Reporting discipline
Honest verified/thin/missing counts per family in every report. Failures reported
as-is; negative results preserved, never re-fished.

## Corpus & compute
MGnify API (16 MAG proteomes: 8 soil-v1-0 + 8 marine-sediment-v1-0, provenance.tsv),
UniProt TrEMBL server-side bg searches (seeded 800-seq subsamples), E. coli K-12 +
B. subtilis decoy proteomes (2-core/2GB-safe). Structure work: colabfold-scale on
selected top candidates only.
