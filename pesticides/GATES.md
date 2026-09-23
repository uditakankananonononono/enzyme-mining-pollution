# LOCKED GATES - Builder 20 slice: pesticide-degrading enzyme mining
# Project 10: enzyme-mining-pollution | Locked 2026-09-23 13:30 IST BEFORE touching outcome data
# Scope: organophosphate (OP) hydrolases + pyrethroid hydrolases from public metagenomes.
# Amended bar (shared w/ builder 19): mining platform + validated specificity gate, NOT another curation DB.

## Reference set (assembled with evidence URLs, locked before mining)
- OP hydrolases across the 3 documented superfamilies (Singh 2009 Nat Rev Microbiol):
  (a) PTE/amidohydrolase: P0A433, P0A434, OpdA (A. radiobacter P230), PteP, plus characterized unreviewed entries;
  (b) MPH family (metallo-beta-lactamase-like): MPH (WBC-3), Mpd (M6), OPHC2;
  (c) PON (P27169 human PON1) and OPAA (Alteromonas Q44238), DFPase.
- Pyrethroid hydrolases from literature with accessions: PytH, PytY, PytZ, Pye3, Sys410,
  Est3385, Est804, EstP, EstSt7, Cest2923, BioH(biobed), mammalian pyrethroid-hydrolyzing carboxylesterases.
- Target: >=30 enzymes per class, each with an experimental-evidence URL (UniProt/BRENDA/paper).
- Split BEFORE mining: seeded 70/30 train/held-out, stratified by family. Split list is frozen.

## GATE 1 - Reference recall >=90% (parent-set)
Pipeline (profile + active-site-motif filter) must recover >=90% of HELD-OUT known enzymes
hidden in a mining corpus (MGnify proteins + UniMES + decoy proteomes).
Held-out sequences are NEVER used to build profiles. Positive-control sanity: 100% recovery
of TRAINING enzymes confirmed unblinded first (pipeline recovers a known answer).
Operating threshold fixed by training-set cross-validation ONLY, before the held-out test. No post-hoc tuning.

## GATE 2 - Specificity discriminator vs documented false-positive modes (parent-set)
Labeled negatives (>=30 per mode):
  (a) non-PTE amidohydrolase-superfamily members incl. phosphotriesterase-like lactonases (PLLs)
      (documented: Afriat et al. 2006; most PTE-like hits are native lactonases);
  (b) esterase-positive but pyrethroid-INACTIVE enzymes from functional screens
      (documented: Sys410 screen 3 est+ clones -> 1 true; Pye3 screen 6 est+ clones -> 1 true);
  (c) generic carboxylesterases/lipases (pNP-ester active, no pesticide evidence).
Required at the frozen operating threshold: precision >=90% on held-out labeled set AND AUROC >=0.9.

## GATE 3 - Novelty
Claimed novel candidates: <80% pairwise identity to any reference enzyme AND absent from
BRENDA/literature curation for OP/pyrethroid activity (checked, not assumed).

## Reporting
Honest verified/thin/missing counts per family in every report. Failures reported as-is; no re-fishing.

## Corpus & compute plan
MGnify sequence-search API + UniProt/UniMES API server-side searches (2-core/2GB sandbox safe).
Structure work: colabfold-scale on selected top candidates only, not hundreds.
